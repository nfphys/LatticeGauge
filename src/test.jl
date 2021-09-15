
function test_initialize_gaugefields()
    param = PhysicalParam()
    @unpack dim, Nsite = param 

    Us_U1 = ones(ComplexF64, Nsite, Nsite, Nsite, Nsite, 4)
    @time initialize_gaugefields!(Us_U1, param; random=true)

    Us_SU2 = Array{SU2}(undef, Nsite, Nsite, Nsite, Nsite, 4)
    @time initialize_gaugefields!(Us_SU2, param; random=true)
end


function test_thermalization_SU2(param; β=1.0, S=1)
    @unpack Nsite, Nthermal, Nsweep = param
    
    #ordered start
    Us_ordered = Array{SU2}(undef, Nsite, Nsite, Nsite, Nsite, 4)
    initialize_gaugefields!(Us_ordered, param)
    
    Ps_ordered = zeros(Float64, Nthermal)
    Ps_average_ordered = zeros(Float64, Nthermal)
    
    Ws_ordered = zeros(Float64, Nthermal)
    Ws_average_ordered = zeros(Float64, Nthermal)
    
    @time for isweep in 1:Nthermal
        heatbath!(Us_ordered, param, β)
        Ps_ordered[isweep] = calc_average_plaquette(param, Us_ordered)
        Ws_ordered[isweep] = calc_wilson_loop(param, Us_ordered, S, S) 
        if isweep === 1
            Ps_average_ordered[isweep] = Ps_ordered[isweep]
            Ws_average_ordered[isweep] = Ws_ordered[isweep]
        else
            Ps_average_ordered[isweep] = Ps_average_ordered[isweep-1]*(isweep-1) + Ps_ordered[isweep]
            Ps_average_ordered[isweep] /= isweep
            
            Ws_average_ordered[isweep] = Ws_average_ordered[isweep-1]*(isweep-1) + Ws_ordered[isweep]
            Ws_average_ordered[isweep] /= isweep
        end
    end
    
    # random start
    Us_random = Array{SU2}(undef, Nsite, Nsite, Nsite, Nsite, 4)
    initialize_gaugefields!(Us_random, param; random=true)
    
    Ps_random = zeros(Float64, Nthermal)
    Ps_average_random = zeros(Float64, Nthermal)
    
    Ws_random = zeros(Float64, Nthermal)
    Ws_average_random = zeros(Float64, Nthermal)
    
    @time for isweep in 1:Nthermal
        heatbath!(Us_random, param, β)
        Ps_random[isweep] = calc_average_plaquette(param, Us_random)
        Ws_random[isweep] = calc_wilson_loop(param, Us_random, S, S) 
        if isweep === 1
            Ps_average_random[isweep] = Ps_random[isweep]
            Ws_average_random[isweep] = Ws_random[isweep]
        else
            Ps_average_random[isweep] = Ps_average_random[isweep-1]*(isweep-1) + Ps_random[isweep]
            Ps_average_random[isweep] /= isweep
            
            Ws_average_random[isweep] = Ws_average_random[isweep-1]*(isweep-1) + Ws_random[isweep]
            Ws_average_random[isweep] /= isweep
        end
    end
    
    # plot 
    p = plot(xlabel="iteration", ylabel="average plaquette", legend=:bottomright)
    plot!(p, Ps_average_ordered; label="ordered start")
    plot!(p, Ps_average_random; label="random start")
    display(p)
    
    p = plot(xlabel="iteration", ylabel="wilson loop")
    plot!(p, Ws_average_ordered; label="ordered start")
    plot!(p, Ws_average_random; label="random start")
    display(p)
end