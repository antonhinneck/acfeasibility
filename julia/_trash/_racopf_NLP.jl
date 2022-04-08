function build_acopf_feas_nlp(_case, pk; Iea = [true for i in 1:length(_case.lines)], x0 = nothing)
    ## Requires Ipopt
    #################

    m = JuMP.Model(with_optimizer(Ipopt.Optimizer))
    # set_optimizer_attribute(m, "ms_enable", 1)
    # set_optimizer_attribute(m, "ms_maxsolves", 30)

    @variable(m, ϵpb[_case.buses] >= 0)
    @variable(m, ϵmb[_case.buses] >= 0)
    # @variable(m, ϵpp[_case.generators] >= 0)
    # @variable(m, ϵmp[_case.generators] >= 0)
    @variable(m, θ[_case.buses])
    @variable(m, γ[l=_case.lines; Iea[l]] >= 0)
    @variable(m, p[_case.generators] >= 0)
    @variable(m, v[_case.buses] >= 0)
    
    # @variable(m, s[l=_case.lines; Iea[l]], Bin)
    #@variable(m, p[g=_case.generators])

    # if x0 != nothing
    #     set_start_value(s[l], x0[l])
    # end

    M = THETAMAX - THETAMIN
    
    #@NLexpression(m, δθ[l = _case.lines; Iea[l]], θ[_case.bus_id[_case.line_end[l]]] - θ[_case.bus_id[_case.line_start[l]]])
    #@expression(m, mδθ[l = _case.lines; Iea[l]], θ[_case.bus_id[_case.line_start[l]]] - θ[_case.bus_id[_case.line_end[l]]])
    @expression(m, bl[l = _case.lines; Iea[l]], 100.0 ./ max(_case.line_x[l], 0.0001))
    @expression(m, gl[l = _case.lines; Iea[l]], 100.0 ./ max(_case.line_r[l], 0.0001))

    @objective(m, Min, sum(ϵpb[b] + ϵmb[b] for b in case.buses))

    #@constraint(m, flowD1[l=_case.lines; Iea[l]], !s[l] => {γ[l] == mδθ[l]})
    #@NLconstraint(m, flowD2[l=_case.lines; Iea[l]], γ[l] == 0.4 * δθ[l]^2)
    @NLexpression(m, vs[l=_case.lines; Iea[l]], v[_case.bus_id[_case.line_end[l]]] * v[_case.bus_id[_case.line_start[l]]])
    @NLconstraint(m, flowP[b = _case.buses], sum( vs[l] * (gl[l] * (1 - 0.4 * (θ[_case.bus_id[_case.line_end[l]]] - θ[_case.bus_id[_case.line_start[l]]])^2) + bl[l] * (θ[_case.bus_id[_case.line_end[l]]] - θ[_case.bus_id[_case.line_start[l]]]))  for l in _case.lines_at_bus[b] if Iea[l])  == _case.bus_Pd[b] / 1.0 + ϵpb[b] - ϵmb[b] - sum( p[g] for g in _case.generators_at_bus[b] ) )

    @constraint(m, θLim1[b = _case.buses], θ[b] <= pi / 3)
    @constraint(m, θLim2[b = _case.buses], θ[b] >= -pi / 3)

    @constraint(m, vLim1[b = _case.buses], v[b] <= _case.bus_Vmax[b])
    @constraint(m, vLim2[b = _case.buses], v[b] >= _case.bus_Vmin[b])

    @constraint(m, genLim1[g = _case.generators], p[g] <= _case.generator_Pmax[g])
    @constraint(m, genLim2[g = _case.generators], p[g] >= _case.generator_Pmin[g])

    # @constraint(m, θLim3[l = _case.lines; Iea[l]], θ[_case.bus_id[_case.line_end[l]]] <= θ[_case.bus_id[_case.line_start[l]]] + s[l] * M)
    # @constraint(m, θLim4[l = _case.lines; Iea[l]], θ[_case.bus_id[_case.line_start[l]]] <= θ[_case.bus_id[_case.line_end[l]]] + (1 - s[l]) * M)

    return m

end