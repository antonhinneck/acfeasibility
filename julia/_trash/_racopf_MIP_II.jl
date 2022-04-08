function build_acopf_feas_mipII(_case,_θ,pk,pi; Iea = [true for i in 1:length(_case.lines)])
    ## Requires Ipopt
    #################

    m = JuMP.Model(with_optimizer(KNITRO.Optimizer))
    #m = JuMP.Model(with_optimizer(AmplNLWriter.Optimizer("/usr/local/bin/bonmin")))
    # m = Model(() -> AmplNLWriter.Optimizer("/usr/local/bin/bonmin"))
    # set_optimizer_attribute(m, "bonmin.nlp_log_level", 1)
    # m = Model(() -> AmplNLWriter.Optimizer("/usr/local/bin/couenne"))
    # set_optimizer_attribute(m, "display_stats", "yes")
    # m = JuMP.Model(with_optimizer(EAGO.Optimizer))
    # set_optimizer_attribute(m, "verbosity", 1)
    #m = JuMP.Model(with_optimizer(Mosek.Optimizer))

    @variable(m, v[_case.buses] >= 0)
    # @variable(m, ϵpp[_case.generators] >= 0)
    # @variable(m, ϵmp[_case.generators] >= 0)
    @variable(m, _p[_case.generators] >= 0)
    @variable(m, ϵpb[_case.buses] >= 0)
    @variable(m, ϵmb[_case.buses] >= 0)

    #@expression(m, sin_δθ[l = _case.lines], sin(_θ[_case.bus_id[_case.line_end[l]]] - _θ[_case.bus_id[_case.line_start[l]]]))
    #@expression(m, cos_δθ[l = _case.lines], cos(_θ[_case.bus_id[_case.line_end[l]]] - _θ[_case.bus_id[_case.line_start[l]]]))
    @expression(m, δθ[l = _case.lines], _θ[_case.bus_id[_case.line_end[l]]] - _θ[_case.bus_id[_case.line_start[l]]])
    @expression(m, bl[l = _case.lines], 1.0 ./ max(_case.line_x[l], 0.001))
    #@expression(m, gl[l = _case.lines], 1 ./ max(_case.line_r[l], 0.001))
    @expression(m, gl[l = _case.lines], 1.0 ./ max(_case.line_r[l], 0.001))
    #@expression(m, vv[l = _case.lines], v[_case.bus_id[_case.line_start[l]]] *  v[_case.bus_id[_case.line_end[l]]] ) 
    # @expression(m, vv[l = _case.lines], v[_case.bus_id[_case.line_start[l]]] *  v[_case.bus_id[_case.line_end[l]]] ) 


    @objective(m, Min, 200 * sum(ϵmb[b] + ϵpb[b] for b in _case.generators) + sum(_p[g] * _case.generator_c1[g] for g in _case.generators))

    @constraint(m, flowP[b = _case.buses], sum( v[_case.bus_id[_case.line_end[l]]] * v[_case.bus_id[_case.line_start[l]]] * (gl[l] * cos(δθ[l]) + bl[l] * sin(δθ[l]))  for l in _case.lines_at_bus[b] if Iea[l])  == _case.bus_Pd[b] + ϵpb[b] - ϵmb[b] - sum( _p[g] for g in _case.generators_at_bus[b] ) )
    #@constraint(m, flowQ[b = _case.buses], ϵmq[b] - ϵpq[b] + sum( vv[l] * (gl[l] * sin_δθ[l] - bl[l] * cos_δθ[l])  for l in _case.lines_at_bus[b] if Iea[l])  == _case.bus_Qd[b] - sum( _q[g] for g in _case.generators_at_bus[b] ) )

    # @constraint(m, reGenLim1[g = _case.generators], _q[g] <= _case.generator_Qmax[g])
    # @constraint(m, reGenLim2[g = _case.generators], _q[g] >= _case.generator_Qmin[g])

    @constraint(m, genLim1[g = _case.generators], _p[g] <= _case.generator_Pmax[g])
    @constraint(m, genLim2[g = _case.generators], _p[g] >= _case.generator_Pmin[g])

    @constraint(m, vLim1[b = _case.buses], v[b] <= _case.bus_Vmax[b])
    @constraint(m, vLim2[b = _case.buses], v[b] >= _case.bus_Vmin[b])

    return m

end