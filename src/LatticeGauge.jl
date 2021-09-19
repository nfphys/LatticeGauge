module LatticeGauge

using Plots
using StaticArrays
using LinearAlgebra
using Parameters 
using Statistics 

include("./U1.jl")
include("./SU2.jl")
include("./SU3.jl")
include("./initialize.jl")
include("./staple.jl")
include("./heatbath.jl")
include("./metropolis.jl")
include("./overrelaxation.jl")
include("./observables.jl")
include("./test.jl")

export measure_observables!

@with_kw struct PhysicalParam 
    dim::Int64 = 4
    Nsite::Int64 = 4 

    Nthermal::Int64 = 100
    Nsweep::Int64 = 100

    do_OR::Bool = false
    NOR::Int64 = 4

    Nhit::Int64 = 10

    method::Symbol = :heatbath 
    @assert method === :heatbath || method === :metropolis
end

function measure_observables!(Us, param, β; random=false)
    @unpack Nsite, Nthermal, Nsweep, do_OR, NOR, method = param 

    initialize_gaugefields!(Us, param; random=random)

    # thermalization
    for isweep in 1:Nthermal
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
    end
    
    # measurement
    Ps = zeros(Float64, Nsweep)
    for isweep in 1:Nthermal
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
        Ps[isweep] = calc_average_plaquette(param, Us)
    end
    
    P  = mean(Ps)
    ΔP = stdm(Ps, P)
    
    return P, ΔP
end




function main!(Us, param, βs)
    @unpack Nsite, dim = param

    @assert ndims(Us) === dim+1
    for i in 1:dim 
        @assert size(Us, i) === Nsite 
    end
    @assert size(Us, dim+1) === dim 

    Ps_ordered = zeros(Float64, length(βs))
    ΔPs_ordered = zeros(Float64, length(βs))

    Ps_random = zeros(Float64, length(βs))
    ΔPs_random = zeros(Float64, length(βs))

    for iβ in 1:length(βs)
        β = βs[iβ]

        # ordered start
        @time P, ΔP = measure_observables!(Us, param, β; random=false)

        Ps_ordered[iβ] = P
        ΔPs_ordered[iβ] = ΔP

        # random start
        @time P, ΔP = measure_observables!(Us, param, β; random=true)
        Ps_random[iβ] = P
        ΔPs_random[iβ] = ΔP
    end

    group = :U1 
    if typeof(Us) == Array{SU2, dim+1} 
        group = :SU2
    end
    if typeof(Us) == Array{SU3, dim+1}
        group = :SU3 
    end

    result = (dim=dim, βs=βs, 
        group=group,
        Ps_ordered=Ps_ordered, 
        ΔPs_ordered=ΔPs_ordered, 
        Ps_random=Ps_random, 
        ΔPs_random=ΔPs_random)

    return result 
end

function plot_result(result; xlog=true, ylog=true, yerr=false, save_figure=false, figure_name="result")
    @unpack dim, βs, group, Ps_ordered, ΔPs_ordered, Ps_random, ΔPs_random = result 

    # plot average plaquette
    p = plot(xlabel="β", ylabel="P", legend=:bottomleft)
    if xlog 
        plot!(p, xaxis=:log)
    end
    if ylog 
        plot!(p, yaxis=:log)
    end
    if yerr
        plot!(p, βs, Ps_ordered; yerr=ΔPs_ordered, marker=:dot, line=false, label="ordered start")
        plot!(p, βs, Ps_random; yerr=ΔPs_random, marker=:dot, line=false, label="random start")
    else
        plot!(p, βs, Ps_ordered; marker=:dot, line=false, label="ordered start")
        plot!(p, βs, Ps_random; marker=:dot, line=false, label="random start")
    end

    function f(x)
        if group === :U1
            return ifelse(x<1.5, 1-x/2 , NaN)
        end
        if group === :SU2 
            return ifelse(x<3.5, 1-x/4 , NaN)
        end
        if group === :SU3 
            return ifelse(x<6.0, 1-x/18-x*x/216, NaN)
        end
        return NaN
    end
    plot!(p, βs, f.(βs), label=false, line=:black)

    function g(x)
        if group === :U1
            return ifelse(x>0.5, 1/(dim*x),  NaN)
        end
        if group === :SU2 
            return ifelse(x>0.8, 3/(dim*x),  NaN)
        end
        if group === :SU3
            return ifelse(x>5.0, 8/(dim*x), NaN)
        end
        return NaN
    end
    plot!(p, βs, g.(βs), label=false, line=:black)
    if save_figure
        savefig("LatticeGauge_figure/average_plaquette_" * figure_name * ".png")
    end
    display(p)


    # plot plaquette fluctuation
    p = plot(xlabel="β", ylabel="ΔP", legend=:bottomleft)
    if xlog 
        plot!(xaxis=:log)
    end
    plot!(p, βs, ΔPs_ordered; marker=:dot, line=false, label="ordered start")
    plot!(p, βs, ΔPs_random; marker=:dot, line=false, label="random start")
    if save_figure
        savefig("LatticeGauge_figure/plaquette_fluctuation_" * figure_name * ".png")
    end
    display(p)
end





end # module
