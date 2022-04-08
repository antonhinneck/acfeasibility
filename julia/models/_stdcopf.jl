function build_stdcopf(ge::Gurobi.Env, data; Iea = [true for i in 1:length(data.lines)], outputflag = 0, threads = 1)

    m = JuMP.Model(with_optimizer(Gurobi.Optimizer, ge))
    set_optimizer_attribute(m, "OutputFlag", outputflag)
    set_optimizer_attribute(m, "Threads", threads)

    @variable(m, p[data.generators] >= 0)
    @variable(m, v[data.buses])
    @variable(m, f[l = data.lines; Iea[l]])

    # JuMP.fix(v[69], 0.0, force = true)

    #Minimal generation costs
    @objective(m, Min, sum(data.generator_c1[g] * p[g] for g in data.generators))

    #Current law
    @constraint(m, nb[n = data.buses],
    sum(p[g] for g in data.generators_at_bus[n]) + sum(f[l] for l in data.lines_start_at_bus[n] if Iea[l]) - sum(f[l] for l in data.lines_end_at_bus[n] if Iea[l]) == data.bus_Pd[n])

    #Voltage law
    @constraint(m, voltage_1[l = data.lines; Iea[l]],
    (data.base_mva / data.line_x[l]) * (v[data.bus_id[data.line_start[l]]] - v[data.bus_id[data.line_end[l]]]) == f[l])

    #Capacity constraint
    @constraint(m, production_capacity1[g = data.generators], p[g] <= data.generator_Pmax[g])
    #@constraint(m, production_capacity2[g = data.generators], p[g] >= 0 )#data.generator_Pmin[g])

    # #Angle limits
    @constraint(m, theta_limit1[n = data.buses], v[n] <= THETAMAX)
    @constraint(m, theta_limit2[n = data.buses], v[n] >= THETAMIN)

    #Line limit
    @constraint(m, fl1[l in data.lines; Iea[l]], f[l] <= data.line_capacity[l])
    @constraint(m, fl2[l in data.lines; Iea[l]], f[l] >= -data.line_capacity[l])
    # @constraint(m, fl3[l in data.lines; !Iea[l]], f[l] == 0.0)

    return m
end
