module LatticeGauge

using Plots
using LinearAlgebra
using Parameters 
using Statistics 

include("./SU2.jl")
include("./initialize.jl")
include("./heatbath.jl")
include("./observables.jl")
include("./test.jl")

export measure_observables!

@with_kw struct PhysicalParam 
    dim::Int64 = 4
    Nsite::Int64 = 4 
    Nthermal::Int64 = 100
    Nsweep::Int64 = 100
    β::Float64 = 1.0
end

function measure_observables!(Us, param, β; random=false)
    @unpack Nsite, Nthermal, Nsweep = param 

    initialize_gaugefields!(Us, param; random=random)

    # thermalization
    for isweep in 1:Nthermal
        heatbath!(Us, param, β)
    end
    
    # measurement
    Ps = zeros(Float64, Nsweep)
    for isweep in 1:Nthermal
        heatbath!(Us, param, β)
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

    group = :SU2 
    if typeof(Us) == Array{ComplexF64, dim+1} 
        group = :U1
    end

    result = (dim=dim, βs=βs, 
        group=group,
        Ps_ordered=Ps_ordered, 
        ΔPs_ordered=ΔPs_ordered, 
        Ps_random=Ps_random, 
        ΔPs_random=ΔPs_random)

    return result 
end

function plot_result(result; xlog=true, yerr=false, save_figure=false, figure_name="result")
    @unpack dim, βs, group, Ps_ordered, ΔPs_ordered, Ps_random, ΔPs_random = result 

    # plot average plaquette
    p = plot(xlabel="β", yaxis=:log, ylabel="P", legend=:bottomleft)
    if xlog 
        plot!(xaxis=:log)
    end
    if yerr
        plot!(p, βs, Ps_ordered; yerr=ΔPs_ordered, marker=:dot, line=false, label="ordered start")
        plot!(p, βs, Ps_random; yerr=ΔPs_random, marker=:dot, line=false, label="random start")
    else
        plot!(p, βs, Ps_ordered; marker=:dot, line=false, label="ordered start")
        plot!(p, βs, Ps_random; marker=:dot, line=false, label="random start")
    end


    function f(x)
        if group === :SU2 
            return ifelse(x<3.5, 1-x/4, NaN)
        end
        if group === :U1
            return ifelse(x<1.5, 1-x/2, NaN)
        end
        return NaN
    end
    plot!(p, βs, f.(βs), label=false, line=:black)

    function g(x)
        if group === :SU2 
            return ifelse(x>0.8, 3/(dim*x),  NaN)
        end
        if group === :U1
            return ifelse(x>0.5, 1/(dim*x),  NaN)
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
