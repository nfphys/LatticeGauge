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


function main_U1(param, βs)
    @unpack Nsite = param
    Us = ones(ComplexF64, Nsite, Nsite, Nsite, Nsite, 4)

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

    result = (βs=βs, 
        Ps_ordered=Ps_ordered, 
        ΔPs_ordered=ΔPs_ordered, 
        Ps_random=Ps_random, 
        ΔPs_random=ΔPs_random)

    return result 
end


function plot_U1(result; xlog=true)
    @unpack βs, Ps_ordered, ΔPs_ordered, Ps_random, ΔPs_random = result 

    # plot average plaquette 
    p = plot(xlabel="β", yaxis=:log, ylabel="P", legend=:bottomleft)
    if xlog 
        plot!(xaxis=:log)
    end
    plot!(p, βs, Ps_ordered; marker=:dot, line=false, label="ordered start")
    plot!(p, βs, Ps_random; marker=:dot, line=false, label="random start")

    f(x) = ifelse(x<1.5, 1-x/2, NaN)
    plot!(p, βs, f.(βs), label="strong coupling limit: 1 - β/2")

    g(x) = ifelse(x>0.5, 1/4x,  NaN)
    plot!(p, βs, g.(βs), label="weak coupling limit: 1/4β")
    display(p)

    # plot fluctuation of plaquette 
    p = plot(xlabel="β", ylabel="ΔP")
    if xlog 
        plot!(xaxis=:log)
    end
    plot!(p, βs, ΔPs_ordered; marker=:dot, line=false, label="ordered start")
    plot!(p, βs, ΔPs_random; marker=:dot, line=false, label="random start")
    display(p)
end





function main_SU2(param, βs)
    @unpack Nsite = param
    Us = Array{SU2}(undef, Nsite, Nsite, Nsite, Nsite, 4)

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

    result = (βs=βs, 
        Ps_ordered=Ps_ordered, 
        ΔPs_ordered=ΔPs_ordered, 
        Ps_random=Ps_random, 
        ΔPs_random=ΔPs_random)

    return result 
end

function plot_SU2(result; xlog=true)
    @unpack βs, Ps_ordered, ΔPs_ordered, Ps_random, ΔPs_random = result 

    # plot average plaquette
    p = plot(xlabel="β", yaxis=:log, ylabel="P", legend=:bottomleft)
    if xlog 
        plot!(xaxis=:log)
    end
    plot!(p, βs, Ps_ordered; marker=:dot, line=false, label="ordered start")
    plot!(p, βs, Ps_random; marker=:dot, line=false, label="random start")


    f(x) = ifelse(x<3.5, 1-x/4, NaN)
    plot!(p, βs, f.(βs), label="strong coupling limit: 1 - β/4")

    g(x) = ifelse(x>0.8, 3/4x,  NaN)
    plot!(p, βs, g.(βs), label="weak coupling limit: 3/4β")
    display(p)

    # plot fluctuation of plaquette
    p = plot(xlabel="β", ylabel="ΔP", legend=:bottomleft)
    if xlog 
        plot!(xaxis=:log)
    end
    plot!(p, βs, ΔPs_ordered; marker=:dot, line=false, label="ordered start")
    plot!(p, βs, ΔPs_random; marker=:dot, line=false, label="random start")
    display(p)
end




function main_SU2_5d(param, βs)
    @unpack Nsite = param
    Us = Array{SU2}(undef, Nsite, Nsite, Nsite, Nsite, Nsite, 5)

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

    result = (βs=βs, 
        Ps_ordered=Ps_ordered, 
        ΔPs_ordered=ΔPs_ordered, 
        Ps_random=Ps_random, 
        ΔPs_random=ΔPs_random)

    return result 
end

function plot_SU2_5d(result; xlog=true)
    @unpack βs, Ps_ordered, ΔPs_ordered, Ps_random, ΔPs_random = result 

    # plot average plaquette
    p = plot(xlabel="β", yaxis=:log, ylabel="P", legend=:bottomleft)
    if xlog 
        plot!(xaxis=:log)
    end
    plot!(p, βs, Ps_ordered; marker=:dot, line=false, label="ordered start")
    plot!(p, βs, Ps_random; marker=:dot, line=false, label="random start")


    f(x) = ifelse(x<3.5, 1-x/4, NaN)
    plot!(p, βs, f.(βs), label="strong coupling limit: 1 - β/4")

    g(x) = ifelse(x>0.8, 3/5x,  NaN)
    plot!(p, βs, g.(βs), label="weak coupling limit: 3/5β")
    display(p)

    # plot fluctuation of plaquette
    p = plot(xlabel="β", ylabel="ΔP", legend=:bottomleft)
    if xlog 
        plot!(xaxis=:log)
    end
    plot!(p, βs, ΔPs_ordered; marker=:dot, line=false, label="ordered start")
    plot!(p, βs, ΔPs_random; marker=:dot, line=false, label="random start")
    display(p)
end






end # module
