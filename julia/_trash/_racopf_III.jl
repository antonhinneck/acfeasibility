function build_acopf_feas_III(_case, _θ; Iea = [true for i in 1:length(_case.lines)])
    ## Requires Ipopt
    #################

    m = JuMP.Model(with_optimizer(KNITRO.Optimizer))
    #m = JuMP.Model(with_optimizer(AmplNLWriter.Optimizer("/usr/local/bin/bonmin")))
    # m = Model(() -> AmplNLWriter.Optimizer("/usr/local/bin/bonmin"))
    # set_optimizer_attribute(m, "bonmin.nlp_log_level", 1)
    # m = Model(() -> AmplNLWriter.Optimizer("/usr/local/bin/couenne"))
    # set_optimizer_attribute(m, "display_stats", "yes")
    # m = JuMP.Model(with_optimizer(EAGO.Optimizer))
    set_optimizer_attribute(m, "ms_enable", 1)
    #m = JuMP.Model(with_optimizer(Mosek.Optimizer))

    @variable(m, v[_case.buses] >= 0)
    @variable(m, em[_case.buses] >= 0)
    @variable(m, ep[_case.buses] >= 0)
    @variable(m, _p[_case.generators] >= 0)
    # @variable(m, _θ[_case.buses] >= 0)

    #@expression(m, sin_δθ[l = _case.lines], sin(_θ[_case.bus_id[_case.line_end[l]]] - _θ[_case.bus_id[_case.line_start[l]]]))
    #@expression(m, cos_δθ[l = _case.lines], cos(_θ[_case.bus_id[_case.line_end[l]]] - _θ[_case.bus_id[_case.line_start[l]]]))
    @expression(m, δθ[l = _case.lines], _θ[_case.bus_id[_case.line_start[l]]] - _θ[_case.bus_id[_case.line_end[l]]])
    @expression(m, bl[l = _case.lines], 100.0 ./ max(_case.line_x[l], 0.001))
    #@expression(m, gl[l = _case.lines], 1 ./ max(_case.line_r[l], 0.001))
    @expression(m, gl[l = _case.lines], 100.0 ./ max(_case.line_r[l], 0.001))
    #@expression(m, vv[l = _case.lines], v[_case.bus_id[_case.line_start[l]]] *  v[_case.bus_id[_case.line_end[l]]] ) 
    # @expression(m, vv[l = _case.lines], v[_case.bus_id[_case.line_start[l]]] *  v[_case.bus_id[_case.line_end[l]]] ) 


    @objective(m, Min, sum(em[b] + ep[b] for b in _case.buses))

    @NLconstraint(m, flowP[b = _case.buses], sum( v[_case.bus_id[_case.line_start[l]]] * v[_case.bus_id[_case.line_end[l]]] * (gl[l] * cos(_θ[_case.bus_id[_case.line_start[l]]] - _θ[_case.bus_id[_case.line_end[l]]]) + bl[l] * sin(_θ[_case.bus_id[_case.line_start[l]]] - _θ[_case.bus_id[_case.line_end[l]]]))  for l in _case.lines_at_bus[b] if Iea[l])  == _case.bus_Pd[b] - sum( _p[g] for g in _case.generators_at_bus[b] ) - ep[b] + em[b] )
    #@constraint(m, flowQ[b = _case.buses], ϵmq[b] - ϵpq[b] + sum( vv[l] * (gl[l] * sin_δθ[l] - bl[l] * cos_δθ[l])  for l in _case.lines_at_bus[b] if Iea[l])  == _case.bus_Qd[b] - sum( _q[g] for g in _case.generators_at_bus[b] ) )

    # @constraint(m, reGenLim1[g = _case.generators], _q[g] <= _case.generator_Qmax[g])
    # @constraint(m, reGenLim2[g = _case.generators], _q[g] >= _case.generator_Qmin[g])

    @constraint(m, reGenLim1[g = _case.generators], _p[g] <= _case.generator_Pmax[g])
    @constraint(m, reGenLim2[g = _case.generators], _p[g] >= _case.generator_Pmin[g])

    # @constraint(m, vLim1[b = _case.buses], _θ[b] <= THETAMAX)
    # @constraint(m, vLim2[b = _case.buses], _θ[b] >= THETAMIN)

    @constraint(m, tLim1[b = _case.buses], v[b] <= _case.bus_Vmax[b])
    @constraint(m, tLim2[b = _case.buses], v[b] >= _case.bus_Vmin[b])

    return m

end