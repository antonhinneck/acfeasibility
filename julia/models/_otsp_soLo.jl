function run_otsp_soLo(data; threads = 4, exporting = true, time_limit = 20, start = false, outputflag = 1, logger_active = false)

    ge = Gurobi.Env()
    num_vars = length(data.generators) + length(data.buses) + length(data.lines) + length(data.lines)

    @inline function restrictiveM()
        out = Dict{Int64, Float64}()
        for i in 1:length(data.lines)
            l = data.lines[i]
            push!(out, l => 2.0 * (data.base_mva / data.line_x[l]) * abs(THETAMAX - THETAMIN) + 1.0)
        end
        return out
    end

    M = restrictiveM()

    m = JuMP.Model(with_optimizer(Gurobi.Optimizer))
    set_optimizer_attribute(m, "TimeLimit", time_limit)
    set_optimizer_attribute(m, "OutputFlag", outputflag)
    set_optimizer_attribute(m, "Threads", threads)

    solution_store = Vector{Vector{Float64}}()

    @variable(m, p[data.generators] >= 0)
    @variable(m, v[data.buses])
    @variable(m, f[data.lines])
    @variable(m, z[data.lines], Bin)

    ## Minimal generation costs
    @objective(m, Min, sum(data.generator_c1[g] * p[g] for g in data.generators))

    ## Current law
    @constraint(m, nb[n = data.buses],
    sum(p[g] for g in data.generators_at_bus[n]) + sum(f[l] for l in data.lines_start_at_bus[n]) - sum(f[l] for l in data.lines_end_at_bus[n]) == data.bus_Pd[n])

    ## Voltage law
    @constraint(m, voltage_1[l = data.lines],
    (data.base_mva / data.line_x[l]) * (v[data.bus_id[data.line_start[l]]] - v[data.bus_id[data.line_end[l]]]) + (1 - z[l]) * M[l] >= f[l])
    
    @constraint(m, voltage_2[l = data.lines],
    (data.base_mva / data.line_x[l]) * (v[data.bus_id[data.line_start[l]]] - v[data.bus_id[data.line_end[l]]]) <= f[l] + (1 - z[l]) * M[l])

    ## Capacity constraint
    @constraint(m, production_capacity1[g = data.generators], p[g] <= data.generator_Pmax[g])
    @constraint(m, production_capacity2[g = data.generators], p[g] >= 0) #data.generator_Pmin[g])

    ## Angle limits
    @constraint(m, theta_limit1[n = data.buses], v[n] <= THETAMAX)
    @constraint(m, theta_limit2[n = data.buses], v[n] >= THETAMIN)

    ## Line limit
    @constraint(m, fl1[l in data.lines], f[l] <=  data.line_capacity[l] * z[l])
    @constraint(m, fl2[l in data.lines], f[l] >= -data.line_capacity[l] * z[l])

    @constraint(m, connected[b in data.buses], sum(z[l] for l in data.lines_at_bus[b]) >= 1)

    @constraint(m, fix[l = data.lines[1:60]], z[l] == 1)
    ## JuMP.fix(v[1], 0.0, force = true)

    @inline function callback_logger(cb_data::CallbackData, cb_where::Cint)

        # Saving progress
        #----------------
        if cb_where == convert(Cint, Gurobi.GRB_CB_MIPNODE)

            if Cint(cb_where) == Gurobi.GRB_CB_MIPNODE
                inc = Ref{Cdouble}()
                lb = Ref{Cdouble}()
                Gurobi.GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_OBJBST, inc)
                Gurobi.GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_OBJBND, lb)
                time = Ref{Cdouble}()
                Gurobi.GRBcbget(cb_data, cb_where, GRB_CB_RUNTIME, time)
                src = 0.0
                save(logger, [inc[], lb[], time[], src])
            end
        end
    end

    @inline function callback_sol_export(cb_data::CallbackData, cb_where::Cint)
        # CALLBACK FUNCTION - HOOKS INTO GUROBI
        #----------------------
        # SOLVER STATES:
        #---------------
            #const CB_POLLING = 0
            #const CB_PRESOLVE = 1
            #const CB_SIMPLEX = 2
            #const CB_MIP = 3
            #const CB_MIPSOL = 4
            #const CB_MIPNODE = 5
            #const CB_MESSAGE = 6
            #const CB_BARRIER = 7

        # Saving progress
        #----------------
        if cb_where == convert(Cint, Gurobi.GRB_CB_MIPSOL)
            # Send current incumbent solution
            #--------------------------------
            solcnt = Ref{Cint}()
            Gurobi.GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_SOLCNT, solcnt) 

            if solcnt[] >= 0
                cinc = Array{Float64}(undef, num_vars)
                Gurobi.GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_SOL, cinc)
                cobj = Ref{Cdouble}()
                Gurobi.GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_OBJBST, cobj)
                cbnd = Ref{Cdouble}()
                Gurobi.GRBcbget(cb_data, cb_where, GRB_CB_MIPSOL_OBJBND, cbnd)
                
                push!(solution_store, [cinc[(switched_idx[1]):(switched_idx[2])]...])
            end
        end
    end

    MOIU.attach_optimizer(m)
    grb_model = backend(m).optimizer.model

    if start
        for i in 1:length(data.lines)
            grb_idx = m.moi_backend.model_to_optimizer_map[index(z[i])].value
            Gurobi.GRBsetdblattrelement(grb_model, "Start", grb_idx - 1, Cdouble(1.0))
        end
    end

    switched_idx = [m.moi_backend.model_to_optimizer_map[index(z[1])].value,
                    m.moi_backend.model_to_optimizer_map[index(z[length(data.lines)])].value]

    user_data = _CallbackUserData(
        grb_model,
        (cb_data, cb_where) -> begin
            callback_sol_export(cb_data, cb_where)
            return
        end,
    )

    if exporting
        ret = GRBsetcallbackfunc(grb_model, grb_callback, user_data)
    end

    GRBupdatemodel(grb_model)

    GRBoptimize(grb_model)

    open(string("logs_sol/118/",data.name,now(),".txt"), "w") do f
        for i in 1:length(solution_store)
            line = ""
            for j in 1:length(solution_store[i])
                if j == length(solution_store[i])
                    line = string(line,round(abs(solution_store[i][j])))
                else
                    line = string(line,round(abs(solution_store[i][j])),",")
                end
            end
            write(f, string(line,"\n"))
        end
    end

    return m
end
