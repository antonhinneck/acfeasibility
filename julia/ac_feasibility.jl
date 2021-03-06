using JuMP, Gurobi, MathOptInterface.Utilities, LinearAlgebra, Mosek, MosekTools, LinearAlgebra
using MathOptInterface, Ipopt ## Cbc, GLPK, HTTP, JSON
#include("C://Users//Anton Hinneck//juliaPackages//GitHub//PowerGrids.jl//src//PowerGrids.jl")
include("/home/anton_hinneck/projects/github/PowerGrids.jl/src/PowerGrids.jl")
using .PowerGrids
using Dates
using Gurobi: CallbackData, _CallbackUserData
using CSV, DataFrames
# using EAGO, AmplNLWriter
# using MosekTools, Mosek
# using KNITRO

cd(@__DIR__)
include("utils/grb_utils.jl")
include("models/_stdcopf.jl")
# include("models/_otsp.jl")
include("models/_otsp_soLo.jl")
# include("models/_racopf_I.jl")
# include("models/_racopf_II.jl")
# include("models/_racopf_III.jl")

#

ge = Gurobi.Env()
THETAMAX = 0.6
THETAMIN = -0.6

PowerGrids.set_csv_path("/home/anton_hinneck/projects/github/pglib2csv/pglib/2020-08-21.19-54-30-275/csv")
PowerGrids.csv_cases(verbose = true)
PowerGrids.select_csv_case(5) # 30as bus
case = PowerGrids.loadCase()

run_otsp_soLo(case, start = true, time_limit = 600)

# logname_sols = "logs_sol/118/pglib_opf_case118_ieee2022-02-19T14:01:41.506.txt"
# data_sols = CSV.read(string(logname_sols), DataFrame, header = 0)

# objectives = Vector{Float64}()

# for i in 1:size(data_sols, 1)
#     Iea = convert.(Bool, [j for j in data_sols[i,:]])
#     m_dcopf = build_stdcopf(ge, case; Iea = Iea, outputflag = 1)
#     optimize!(m_dcopf)
#     push!(objectives, objective_value(m_dcopf))
# end

# open("1354_dc_feasibility.csv", "w") do io
#     line = ""
#     ctr = 1
#     for i in objectives
#         line = string(line, i)
#         if ctr < length(objectives)
#             line = string(line, ",")
#         end
#         ctr += 1
#     end
#     write(io, string(line, "\n"))
# end;
















#################################################
sols = DataFrame(CSV.File(string("logs_sol/118/pglib_opf_case118_ieee2022-02-19T14:00:35.049.txt"), skipto = 1))

objs = Vector{Float64}()
vio_p = Vector{Vector{Float64}}()
vio_m = Vector{Vector{Float64}}()

#Iea = convert.(Bool,[sols[11,:]...])
Iea = [true for i in 1:length(case.lines)]
m_dcopf = build_stdcopf(ge, case; Iea = Iea, outputflag = 1)
optimize!(m_dcopf)
pk = value.(m_dcopf[:p]).data
??k = value.(m_dcopf[:v]).data
obj = 0
for i in case.generators
    obj += case.generator_c1[i] * pk[i]
end



include("models/_racopf_NLP.jl")
m_ac_feas_nlp = build_acopf_feas_nlp(case, ??k; Iea)
optimize!(m_ac_feas_nlp)
ipoptmodel = backend(m_ac_feas_nlp).optimizer.model.inner
m_ac_feas_nlp.moi_backend.model_to_optimizer_map

using Ipopt, JuMP

#m = JuMP.Model(with_optimizer(KNITRO.Optimizer))
m = JuMP.Model(with_optimizer(Ipopt.Optimizer))
# set_optimizer_attribute(m, "print_level", 5)
# set_optimizer_attribute(m, "tol", 1e-6)
# set_optimizer_attribute(m, "max_iter", 2400)
# set_optimizer_attribute(m, "dual_inf_tol", 1e-6)
# set_optimizer_attribute(m, "compl_inf_tol", 1e-6)
# set_optimizer_attribute(m, "acceptable_tol", 1e-6)
# set_optimizer_attribute(m, "acceptable_compl_inf_tol", 1e-6)
# set_optimizer_attribute(m, "mu_strategy", "adaptive")

@variable(m, ??pb[case.buses] >= 0)
@variable(m, ??mb[case.buses] >= 0)
@variable(m, ??[case.buses])
#@variable(m, ??[l=case.lines; Iea[l]] >= 0)
@variable(m, p[case.generators] >= 0)
@variable(m, v[case.buses] >= 0)

@expression(m, bl[l = case.lines; Iea[l]], case.base_mva ./ max(case.line_x[l], 0.0001))
@expression(m, gl[l = case.lines; Iea[l]], case.base_mva ./ max(case.line_r[l], 0.0001))

@NLexpression(m, ????[l=case.lines; Iea[l]], ??[case.bus_id[case.line_end[l]]] - ??[case.bus_id[case.line_start[l]]])
@NLexpression(m, vs1[l=case.lines; Iea[l]], v[case.bus_id[case.line_end[l]]] * v[case.bus_id[case.line_start[l]]])
@NLexpression(m, vs2[l=case.lines; Iea[l]], v[case.bus_id[case.line_start[l]]] * v[case.bus_id[case.line_start[l]]])
@NLexpression(m, pi[l=case.lines; Iea[l]], vs[l] * ( gl[l] * cos(????[l]) + bl[l] * sin(????[l])) + gl[l] * vs2[l] )
#@expression(m, m????[l = _case.lines; Iea[l]], ??[_case.bus_id[_case.line_start[l]]] - ??[_case.bus_id[_case.line_end[l]]])
# @expression(m, bl[l = _case.lines; Iea[l]], 1.0 ./ max(_case.line_x[l], 0.0001))
# @expression(m, gl[l = _case.lines; Iea[l]], 1.0 ./ max(_case.line_r[l], 0.0001))

@objective(m, Min, sum(??pb[b] + ??mb[b] for b in case.buses))

#@constraint(m, flowD1[l=_case.lines; Iea[l]], !s[l] => {??[l] == m????[l]})
##@NLconstraint(m, flowD2[l=_case.lines; Iea[l]], ??[l] == 0.4 * ????[l]^2)

#@NLconstraint(m, flowP[b=case.buses], sum( vs[l] * (gl[l] * (1 - 0.4 * (??[case.bus_id[case.line_end[l]]] - ??[case.bus_id[case.line_start[l]]])^2) + bl[l] * (??[case.bus_id[case.line_end[l]]] - ??[case.bus_id[case.line_start[l]]]))  for l in case.lines_at_bus[b] if Iea[l])  == case.bus_Pd[b] + ??pb[b] - ??mb[b] - sum( p[g] for g in case.generators_at_bus[b] ) )
@NLconstraint(m, balance[b=case.buses], sum( pi[l] for l in case.lines_end_at_bus[case.bus_id[b]] if Iea[l]) - sum( pi[l] for l in case.lines_start_at_bus[case.bus_id[b]] if Iea[l])  == sum( p[g] for g in case.generators_at_bus[b] ) - case.bus_Pd[b] / case.base_mva )

@constraint(m, ??Lim1[b = case.buses], ??[b] <= 1.2)
@constraint(m, ??Lim2[b = case.buses], ??[b] >= -1.2)

@constraint(m, vLim1[b = case.buses], v[b] <= case.bus_Vmax[b])
@constraint(m, vLim2[b = case.buses], v[b] >= case.bus_Vmin[b])

#@constraint(m, objc, sum(??mb[b] + ??pb[b] for b in case.buses) <= 0.4)
#@constraint(m, objc[b = case.buses], ??mb[b] + ??pb[b] <= 0.04)

@constraint(m, genLim1[g = case.generators], p[g] <= case.generator_Pmax[g] / case.base_mva)
@constraint(m, genLim2[g = case.generators], p[g] >= case.generator_Pmin[g] / case.base_mva)

# Ipopt.setwarmstart!(ipoptmodel, x0)
# Ipopt.getobjval(ipoptmodel)
# Ipopt.IpoptSolve(ipoptmodel)
# Ipopt.IpoptMathProgModel

#x0 = [0.01 for i in 1:673]
optimize!(m)

ipoptmodel = backend(m).optimizer.model.inner
sol = ipoptmodel.x
x0 = 0.1 .* sol
# x0 = [zeros(length(case.buses))...,
#       zeros(length(case.buses))...,
#       0.1 .+ zeros(length(case.buses))...,
#       #(??k)...,
#     #   [1 for i in 1:length(case.lines) if Iea[i]]...,
#       pk...,
#       #zeros(length(case.generators))...,
#       #[case.generator_Pmax[g] * 0.00001 for g in case.generators]...,
#       [0.96 for i in 1:length(case.buses)]...
# ]
# ipoptmodel.x = x0
#ipoptmodel.x = [0.1 for i in 1:673]
Ipopt.solveProblem(ipoptmodel)

# @variable(m, ??pb[case.buses] >= 0)
# @variable(m, ??mb[case.buses] >= 0)
# @variable(m, ??[case.buses])
# @variable(m, p[case.generators] >= 0)
# @variable(m, v[case.buses] >= 0)

# for i in 1:size(sols, 1)
#     Iea = convert.(Bool,[sols[i,:]...])
#     m_dcopf = build_stdcopf(ge, case; Iea = Iea, outputflag = 1)

#     include("models/_racopf_NLP.jl")
#     m_ac_feas_nlp = build_acopf_feas_nlp(case, pk; Iea)
#     optimize!(m_ac_feas_nlp)

#     push!(vio_p, value.(m_ac_feas_nlp[:??pb]).data)
#     push!(vio_m, value.(m_ac_feas_nlp[:??mb]).data)
#     push!(objs, objective_value(m_ac_feas_nlp))
# end

include("models/_racopf_NLP.jl")
m_ac_feas_nlp = build_acopf_feas_nlp(case, pk; Iea)
optimize!(m_ac_feas_nlp)

diff = Vector{Float64}()
for g in case.generators
    push!(diff, case.generator_Pmax[g] - pk[g])
end
minimum(diff)

print(value.(m_ac_feas_nlp[:gl]).data - value.(m_ac_feas_nlp[:bl]).data)
print(value.(m_ac_feas_nlp[:??]).data)
approx = Vector{Float64}()
exact = Vector{Float64}()
diff  = Vector{Float64}()
for l in case.lines
    if Iea[l]
        push!(approx, 1 - value(m_ac_feas_nlp[:??][l]))
        push!(exact, cos(value(m_ac_feas_nlp[:????][l])))
        push!(diff, value(m_ac_feas_nlp[:gl][l]) - value(m_ac_feas_nlp[:bl][l]))
    end
end

exact
approx
diff

include("models/_racopf_MIP.jl")
m_ac_feas_mip = build_acopf_feas_mip(case, pk; Iea)
optimize!(m_ac_feas_mip)
print(value.(m_ac_feas_mip[:??]).data)
print(value.(m_ac_feas_mip[:??]).data)
print(value.(m_ac_feas_mip[:s]).data)
print(value.(m_ac_feas_mip[:??mp]).data)
print(value.(m_ac_feas_mip[:??pp]).data)

value(m_ac_feas_mip[:m????][1])
value(m_ac_feas_mip[:????][1])
value(m_ac_feas_mip[:s][1])
# value(m_ac_feas_mip[:m????])

vals_approx = Vector{Float64}()
vals_exact = Vector{Float64}()
vals = Vector{Float64}()
for l in case.lines
    if Iea[l]
        push!(vals_exact, cos(value(m_ac_feas_mip[:????][l])))
        if value(m_ac_feas_mip[:s][l]) == 0.0
            push!(vals_approx, 1 - value(m_ac_feas_mip[:m????][l]))
        else
            push!(vals_approx, 1 - value(m_ac_feas_mip[:????][l]))
        end
        push!(vals, value(m_ac_feas_mip[:????][l]))
    end
end

vals_approx
vals_exact
vals

sum(value.(m_ac_feas_mip[:??pp]).data)
sum(value.(m_ac_feas_mip[:??mp]).data)
pk += (value.(m_ac_feas_mip[:??pp]).data - value.(m_ac_feas_mip[:??mp]).data)
??k = value.(m_ac_feas_mip[:??]).data
bi = (value.(m_ac_feas_mip[:??pb]).data - value.(m_ac_feas_mip[:??mb]).data)
obj = 0.0
for i in case.generators
    obj += case.generator_c1[i] * pk[i]
end
obj

include("models/_racopf_MIP_II.jl")
m_ac_feas_mipII = build_acopf_feas_mipII(case, ??k, pk, bi; Iea)
optimize!(m_ac_feas_mipII)
sum(value.(m_ac_feas_mipII[:??pb]).data)
sum(value.(m_ac_feas_mipII[:??mb]).data)
sum([values(case.bus_Pd)...])

include("models/_racopf_I.jl")
m_ac_feas_I = build_acopf_feas_I(case, pk; Iea)
optimize!(m_ac_feas_I)
print(value.(m_ac_feas_I[:??]).data)
println(objective_value(m_ac_feas_I))
sum(value.(m_ac_feas_I[:??pp]).data)
sum(value.(m_ac_feas_I[:??mp]).data)

println(m_ac_feas_I[:flowP][1])
println(round.(value.(m_ac_feas_I[:??]).data, digits = 4))

??k = value.(m_ac_feas_I[:??]).data

value(m_ac_feas_I[:bl][1]) * value(m_ac_feas_I[:????][1])
value(m_ac_feas_I[:????][3])
value(m_ac_feas_I[:abs_c????][3])
value(m_ac_feas_I[:gl][1]) * (1 - value(m_ac_feas_I[:abs_c????][1])) + value(m_ac_feas_I[:bl][1]) * value(m_ac_feas_I[:????][1])

a = [value(m_ac_feas_I[:abs_c????][l]) for l in case.lines if Iea[l]]
b = [value(m_ac_feas_I[:????][l]) for l in case.lines if Iea[l]]
round.(a.- abs.(b), digits = 4)

include("models/_racopf_II.jl")
m_ac_feas_II = build_acopf_feas_II(case, pk, ??k; Iea)
println(m_ac_feas_II[:flowP][1])
optimize!(m_ac_feas_II)
println(m_ac_feas_I[:flowP])
println(round.(value.(m_ac_feas_I[:??]).data, digits = 4))
println(objective_value(m_ac_feas_I))
println(termination_status(m_ac_feas_I))

include("models/_racopf_III.jl")
m_ac_feas_III = build_acopf_feas_III(case, ??k; Iea)
optimize!(m_ac_feas_III)

println(value.(m_ac_feas_II[:????]).data[1:5])
println(value.(m_ac_feas_II[:abs_c????]).data[1:5])
# println(value.(m_ac_feas_II[:abs_c????]).data)
# println(sum(value.(m_ac_feas_II[:abs_c????]).data))
println(sum(value.(m_ac_feas_II[:p]).data))
#println(sum(value.(m_ac_feas_II[:pp]).data))
#println(round.(value.(m_ac_feas_II[:pm]).data))
println(round.(value.(m_ac_feas_II[:p]).data))
# println(sum(value.(m_ac_feas_II[:??mq]).data))
# println(sum(value.(m_ac_feas_II[:??pq]).data))
println(objective_value(m_ac_feas_II))

??k = value.(m_ac_feas_II[:??]).data
qk = value.(m_ac_feas_II[:q]).data
println(??k)
println(qk)

m_ac_feas_I = build_acopf_feas_I(case,p,??k,qk; Iea)
print(m_ac_feas_I)
optimize!(m_ac_feas_I)
objective_value(m_ac_feas_I)