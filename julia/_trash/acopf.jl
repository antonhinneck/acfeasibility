using JuMP, Gurobi, MathOptInterface.Utilities, LinearAlgebra, Mosek, MosekTools, LinearAlgebra
using MathOptInterface, Ipopt ## Cbc, GLPK, HTTP, JSON
#include("C://Users//Anton Hinneck//juliaPackages//GitHub//PowerGrids.jl//src//PowerGrids.jl")
include("/home/anton_hinneck/projects/github/PowerGrids.jl/src/PowerGrids.jl")
using .PowerGrids
using Dates
using Gurobi: CallbackData, _CallbackUserData
using CSV, DataFrames
using EAGO, AmplNLWriter
using MosekTools, Mosek
using KNITRO

cd(@__DIR__)
include("utils/grb_utils.jl")
include("models/_stdcopf.jl")
include("models/_otsp.jl")
include("models/_otsp_soLo.jl")
include("models/_racopf_I.jl")
include("models/_racopf_II.jl")
include("models/_racopf_III.jl")

#run_otsp_soLo(case)

ge = Gurobi.Env()
THETAMAX = 0.6
THETAMIN = -0.6

PowerGrids.set_csv_path("/home/anton_hinneck/projects/github/pglib2csv/pglib/2020-08-21.19-54-30-275/csv")
PowerGrids.csv_cases(verbose = true)
PowerGrids.select_csv_case(4) # 30as bus
case = PowerGrids.loadCase()

Iea = [true for i in case.lines]
sols = DataFrame(CSV.File(string("logs_sol/118/pglib_opf_case118_ieee2022-02-19T14:00:35.049.txt"), skipto = 1))


m = JuMP.Model(with_optimizer(Ipopt.Optimizer))
# set_optimizer_attribute(m, "ms_enable", 1)
# set_optimizer_attribute(m, "ms_maxsolves", 4)
# m = JuMP.Model(with_optimizer(Ipopt.Optimizer))
# set_optimizer_attribute(m, "tol", 1e-1)
# set_optimizer_attribute(m, "constr_viol_tol", 1e-1)
# set_optimizer_attribute(m, "acceptable_constr_viol_tol", 1e-1)
# set_optimizer_attribute(m, "dual_inf_tol", 1e-1)
# set_optimizer_attribute(m, "compl_inf_tol", 1e-1)
# set_optimizer_attribute(m, "mu_init", 1e-3)
# set_optimizer_attribute(m, "s_max", 1e+1)
Iea = convert.(Bool,[sols[7,:]...])
m = JuMP.Model(with_optimizer(Ipopt.Optimizer))
set_optimizer_attribute(m, "acceptable_constr_viol_tol", 1e+0)
@variable(m, p[case.generators] >= 0)
# @variable(m, f_p[l=case.lines; Iea[l]])
# @variable(m, θ[case.buses])
@variable(m, vj[case.buses] >= 0)
@variable(m, vr[case.buses] >= 0)

@objective(m, Min, sum(p[g] * case.generator_c1[g] for g in case.generators) )

@expression(m, bl[l = case.lines; Iea[l]], 1 / max(case.line_x[l], 0.5))
@expression(m, gl[l = case.lines; Iea[l]], 1 / max(case.line_r[l], 0.5))
# @expression(m, gs[b=case.buses], sum(gl[l] for l in case.lines_start_at_bus[b]) - sum(gl[l] for l in case.lines_end_at_bus[b]))
# @NLexpression(m, δθ[l=case.lines; Iea[l]], θ[case.bus_id[case.line_start[l]]] - θ[case.bus_id[case.line_end[l]]])
@NLexpression(m, sum_g[b=case.buses], sum(p[g] for g in case.generators_at_bus[b]) )
# @NLexpression(m, sum_f_p[b=case.buses], sum(f_p[b] for l in case.lines_at_bus[b]) )
# @NLexpression(m, vs1[l=case.lines; Iea[l]], v[case.bus_id[case.line_end[l]]] * v[case.bus_id[case.line_start[l]]])
# @NLexpression(m, vs2[l=case.lines; Iea[l]], v[case.bus_id[case.line_start[l]]] * v[case.bus_id[case.line_start[l]]])
@expression(m, lf[l=case.lines; Iea[l]], case.bus_id[case.line_start[l]])
@expression(m, lt[l=case.lines; Iea[l]], case.bus_id[case.line_end[l]])
@NLexpression(m, rect_1[b=case.buses], sum(gl[l] * vj[lt[l]] - bl[l] * vr[lt[l]] for l in case.lines_at_bus[b] if Iea[l]) )
@NLexpression(m, rect_2[b=case.buses], sum(gl[l] * vr[lt[l]] + bl[l] * vj[lt[l]] for l in case.lines_at_bus[b] if Iea[l]) )

# @NLconstraint(m, f_p_1[l=case.lines; Iea[l]], f_p[l] == vs1[l] * (gl[l] * cos(δθ[l]) * vr[lt[l]]  + bl[l] * sin(δθ[l]))) + bl[l] * vj[lt[l]] + bl[l] * vr[lt[l]]
# @NLconstraint(m, f_p_1[l=case.lines; Iea[l]], f_p[l] == vr[lf[l]] * rect_1[l] + vj[lf[l]] * rect_2[l])
@NLconstraint(m, nb[b=case.buses], sum_g[b] - case.bus_Pd[b] == vr[b] * (rect_1[b]) + vj[b] * (rect_2[b])) 
# @NLconstraint(m)

# @NLconstraint(m, v1[b=case.buses], case.bus_Vmin[b]^2 <= vr[b]^2 + vj[b]^2)
# @NLconstraint(m, v2[b=case.buses], vr[b]^2 + vr[b]^2 <= case.bus_Vmax[b]^2)
# @NLconstraint(m, v3[b=case.buses], case.bus_Vmin[b]^2 <= vj[b]^2 + vj[b]^2)
# @NLconstraint(m, v4[b=case.buses], vr[b]^2 + vj[b]^2 <= case.bus_Vmax[b]^2)
@NLconstraint(m, v1[b=case.buses], 0.94 <= vr[b]^2 + vj[b]^2)
@NLconstraint(m, v2[b=case.buses], vr[b]^2 + vj[b]^2 <= 1.06)

optimize!(m)

obj = 0
for g in case.generators
    obj += value(m[:p][g]) * case.generator_c1[g]
end
println(obj)

ipoptmodel = backend(m).optimizer.model.inner
ipoptmodel.x = [zeros(length(case.generators))..., ones(length(case.buses))..., ones(length(case.buses))...]
Ipopt.solveProblem(ipoptmodel)
sol = backend(m).optimizer.model.inner.x
#Ipopt.solveProblem(ipoptmodel)
# @constraint(m, θLim1[b = case.buses], θ[b] <= 1.04)
# @constraint(m, θLim2[b = case.buses], θ[b] >= -1.04)

# @constraint(m, vLim1[b = case.buses], v[b] <= case.bus_Vmax[b])
# @constraint(m, vLim2[b = case.buses], v[b] >= case.bus_Vmin[b])

# @constraint(m, pLim1[g = case.generators], p[g] <= case.generator_Pmax[g] / case.base_mva)
# @constraint(m, pLim2[g = case.generators], p[g] >= case.generator_Pmin[g] / case.base_mva)

