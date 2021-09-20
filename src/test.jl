
function test_initialize_gaugefields()
    param = PhysicalParam()
    @unpack dim, Nsite = param 

    Us_U1 = ones(ComplexF64, Nsite, Nsite, Nsite, Nsite, 4)
    @time initialize_gaugefields!(Us_U1, param; random=true)

    Us_SU2 = Array{SU2}(undef, Nsite, Nsite, Nsite, Nsite, 4)
    @time initialize_gaugefields!(Us_SU2, param; random=true)

    Us_SU3 = Array{SU3}(undef, Nsite, Nsite, Nsite, Nsite, 4)
    @time initialize_gaugefields!(Us_SU3, param; random=true)
end

function test_calc_average_plaquette()
    param = PhysicalParam()
    @unpack dim, Nsite = param 

    Us_U1 = ones(ComplexF64, Nsite, Nsite, Nsite, Nsite, 4)
    @time initialize_gaugefields!(Us_U1, param; random=true)
    @time P_U1 = calc_average_plaquette(param, Us_U1)
    @show P_U1

    Us_SU2 = Array{SU2}(undef, Nsite, Nsite, Nsite, Nsite, 4)
    @time initialize_gaugefields!(Us_SU2, param; random=true)
    @time P_SU2 = calc_average_plaquette(param, Us_SU2)
    @show P_SU2

    Us_SU3 = Array{SU3}(undef, Nsite, Nsite, Nsite, Nsite, 4)
    @time initialize_gaugefields!(Us_SU3, param; random=true)
    @time P_SU3 = calc_average_plaquette(param, Us_SU3)
    @show P_SU3

    return 
end



function test_thermalization!(Us, param; β=1.0, Nthermal=100)
    @unpack Nsite, method, do_OR, NOR = param
    
    #ordered start
    initialize_gaugefields!(Us, param)
    
    Ps_ordered = zeros(Float64, Nthermal)
    Ps_average_ordered = zeros(Float64, Nthermal)
    Ps_ordered[1] = calc_average_plaquette(param, Us)
    Ps_average_ordered[1] = Ps_ordered[1]
    
    @time for isweep in 2:Nthermal
        if method === :heatbath
            heatbath!(Us, param, β)
        end 
        if method === :metropolis 
            metropolis!(Us, param, β)
        end

        if do_OR
            for iOR in 1:NOR
                overrelaxation!(Us, param, β)
            end
        end

        Ps_ordered[isweep] = calc_average_plaquette(param, Us)
        Ps_average_ordered[isweep] = Ps_average_ordered[isweep-1]*(isweep-1) + Ps_ordered[isweep]
        Ps_average_ordered[isweep] /= isweep
    end
    
    # random start
    initialize_gaugefields!(Us, param; random=true)
    
    Ps_random = zeros(Float64, Nthermal)
    Ps_average_random = zeros(Float64, Nthermal)
    Ps_random[1] = calc_average_plaquette(param, Us)
    Ps_average_random[1] = Ps_random[1]
    
    @time for isweep in 2:Nthermal
        if method === :heatbath
            heatbath!(Us, param, β)
        end 
        if method === :metropolis 
            metropolis!(Us, param, β)
        end

        if do_OR
            for iOR in 1:NOR
                overrelaxation!(Us, param, β)
            end
        end

        Ps_random[isweep] = calc_average_plaquette(param, Us)
        Ps_average_random[isweep] = Ps_average_random[isweep-1]*(isweep-1) + Ps_random[isweep]
        Ps_average_random[isweep] /= isweep
    end
    
    # plot 
    p = plot(xlabel="iteration", ylabel="average plaquette", legend=:bottomright)
    plot!(p, Ps_average_ordered; label="ordered start")
    plot!(p, Ps_average_random; label="random start")
    display(p)
end