function build_acopf_feas_I(_case, pk; Iea = [true for i in 1:length(_case.lines)])
    ## Requires Ipopt
    #################

    m = JuMP.Model(with_optimizer(KNITRO.Optimizer))
    #set_optimizer_attribute(m, "ms_enable", 1)

    #@variable(m, p[_case.generators] >= 0)
    @variable(m, ϵpp[_case.generators] >= 0)
    @variable(m, ϵmp[_case.generators] >= 0)
    @variable(m, abs_cδθ[l=_case.lines; Iea[l]] >= 0)
    @variable(m, θ[_case.buses] >= 0)
    # @variable(m, v[l=_case.lines; Iea[l]] >= 0)
    #@variable(m, q[_case.buses])

    #@expression(m, sin_δθ[l = _case.lines], sin(θ[_case.bus_id[_case.line_end[l]]] - θ[_case.bus_id[_case.line_start[l]]]))
    #@expression(m, cos_δθ[l = _case.lines], cos(θ[_case.bus_id[_case.line_end[l]]] - θ[_case.bus_id[_case.line_start[l]]]))
    # @expression(m, sδθ[l = _case.lines], θ[_case.bus_id[_case.line_start[l]]] - θ[_case.bus_id[_case.line_end[l]]])
    # @expression(m, cδθ[l = _case.lines], θ[_case.bus_id[_case.line_start[l]]] - θ[_case.bus_id[_case.line_end[l]]])
    @NLexpression(m, δθ[l = _case.lines], θ[_case.bus_id[_case.line_end[l]]] - θ[_case.bus_id[_case.line_start[l]]])
    @expression(m, bl[l = _case.lines], 1.0 ./ max(_case.line_x[l], 0.001))
    #@expression(m, gl[l = _case.lines], 1 ./ max(_case.line_r[l], 0.001))
    @expression(m, gl[l = _case.lines], 1.0 ./ max(_case.line_r[l], 0.001))
    #@expression(m, vv[l = _case.lines], v[_case.bus_id[_case.line_start[l]]] *  v[_case.bus_id[_case.line_end[l]]] ) 
    #@expression(m, vv[l = _case.lines], 1.0 )

    # + 500 * sum(abs_cδθ[l] for l in _case.lines if Iea[l])
    @objective(m, Min, sum(abs_cδθ[l] for l in _case.lines if Iea[l]) + sum(ϵpp[g] + ϵmp[g] for g in case.generators))
    #@constraint(m, sum(abs_cδθ[l] for l in _case.lines if Iea[l]) <= 126.0)
    #@objective(m, Min, sum(ϵpp[b] + ϵmp[b] + ϵpq[b] + ϵmq[b] for b in _case.buses) + sum(abs_cδθ[l] for l in _case.lines))
    @NLconstraint(m, abs1[l = _case.lines; Iea[l]], abs_cδθ[l] >= δθ[l])
    @NLconstraint(m, abs2[l = _case.lines; Iea[l]], abs_cδθ[l] >= -δθ[l])
    @NLconstraint(m, flowP[b = _case.buses],  sum( (gl[l] * sqrt(δθ[l] * δθ[l]) + bl[l] * δθ[l])  for l in _case.lines_at_bus[b] if Iea[l])  == _case.bus_Pd[b] / 100.0 - sum( pk[g] + ϵpp[g] - ϵmp[g] for g in _case.generators_at_bus[b] ) )
    #@constraint(m, flowP[b = _case.buses],  sum( (gl[l] * cos(δθ[l]) + bl[l] * sin(δθ[l]))  for l in _case.lines_at_bus[b] if Iea[l])  == _case.bus_Pd[b] - sum( pk[g] for g in _case.generators_at_bus[b] ) )
    #@constraint(m, flowP[b = _case.buses], ϵpp[b] - ϵmp[b] + sum( (gl[l] * abs_cδθ[l] + bl[l] * sδθ[l])  for l in _case.lines_at_bus[b] if Iea[l])  == _case.bus_Pd[b] - sum( p[g] for g in _case.generators_at_bus[b] ) )
    #@constraint(m, flowQ[b = _case.buses], ϵpq[b] - ϵmq[b] + sum( (gl[l] * sδθ[l] - bl[l] * abs_cδθ[l])  for l in _case.lines_at_bus[b] if Iea[l])  == _case.bus_Qd[b] - sum( q[g] for g in _case.generators_at_bus[b] ) )

    # @constraint(m, reGenLim1[g = _case.generators], q[g] <= _case.generator_Qmax[g])
    # @constraint(m, reGenLim2[g = _case.generators], q[g] >= _case.generator_Qmin[g])

    # @constraint(m, genLim1[g = _case.generators], pk[g] + ϵpp[g] <= _case.generator_Pmax[g])
    # @constraint(m, genLim2[g = _case.generators], pk[g] - ϵmp[g] >= _case.generator_Pmin[g])

    @constraint(m, θLim1[b = _case.buses], θ[b] <= 0.1)
    @constraint(m, θLim2[b = _case.buses], θ[b] >= -0.1)
    # @constraint(m, vLim1[l = _case.lines; Iea[l]], v[l] <= _case.bus_Vmax[1]^2)
    # @constraint(m, vLim2[l = _case.lines; Iea[l]], v[l] >= _case.bus_Vmin[1]^2)

    return m

end