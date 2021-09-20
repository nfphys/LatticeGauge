"""
    calc_plaquette(Us, n₁, n₂, μ, ν)

Calculate an ordered-product along the plaquette (n, μ, ν) in 2 dimension.
"""
function calc_plaquette(Us, n₁, n₂, μ, ν)
    Nsite = size(Us, 1)

    δ(μ, ν) = kronecker_delta(μ, ν)
    
    f(n) = check_boundary(n, Nsite)

    U₁ = Us[n₁, n₂, μ]
    U₂ = Us[f(n₁+δ(1,μ)), f(n₂+δ(2,μ)), ν]
    U₃ = conj(Us[f(n₁+δ(1,ν)), f(n₂+δ(2,ν)), μ])
    U₄ = conj(Us[n₁, n₂, ν])

    U = U₁*U₂*U₃*U₄
end


"""
    calc_plaquette(Us, n₁, n₂, n₃, μ, ν)

Calculate an ordered-product along the plaquette (n, μ, ν) in 3 dimension.
"""
function calc_plaquette(Us, n₁, n₂, n₃, μ, ν)
    Nsite = size(Us, 1)

    δ(μ, ν) = kronecker_delta(μ, ν)
    
    f(n) = check_boundary(n, Nsite)

    U₁ = Us[n₁, n₂, n₃, μ]
    U₂ = Us[f(n₁+δ(1,μ)), f(n₂+δ(2,μ)), f(n₃+δ(3,μ)), ν]
    U₃ = conj(Us[f(n₁+δ(1,ν)), f(n₂+δ(2,ν)), f(n₃+δ(3,ν)), μ])
    U₄ = conj(Us[n₁, n₂, n₃, ν])

    U = U₁*U₂*U₃*U₄
end


"""
    calc_plaquette(Us, n₁, n₂, n₃, n₄, μ, ν)

Calculate an ordered-product along the plaquette (n, μ, ν) in 4 dimension.
"""
function calc_plaquette(Us, n₁, n₂, n₃, n₄, μ, ν)
    Nsite = size(Us, 1)

    δ(μ, ν) = kronecker_delta(μ, ν)
    
    f(n) = check_boundary(n, Nsite)

    U₁ = Us[n₁, n₂, n₃, n₄, μ]
    U₂ = Us[f(n₁+δ(1,μ)), f(n₂+δ(2,μ)), f(n₃+δ(3,μ)), f(n₄+δ(4,μ)), ν]
    U₃ = conj(Us[f(n₁+δ(1,ν)), f(n₂+δ(2,ν)), f(n₃+δ(3,ν)), f(n₄+δ(4,ν)), μ])
    U₄ = conj(Us[n₁, n₂, n₃, n₄, ν])

    U = U₁*U₂*U₃*U₄
end



"""
    calc_average_plaquette(param, Us::Array{T, 3}) where T

Calculate average plaquette for gauge fields in 2 dimension.
"""
function calc_average_plaquette(param, Us::Array{T, 3}) where T
    @unpack Nsite = param

    Nc = 1
    if T === SU2 
        Nc = 2
    end
    if T === SU3 
        Nc = 3
    end
    
    P = 0.0
    for ν in 1:2, μ in 1:ν-1, n₂ in 1:Nsite, n₁ in 1:Nsite
        U = calc_plaquette(Us, n₁, n₂, μ, ν)
        P += 1 - real(tr(U))/Nc
    end
    P /= Nsite^2
    return P
end


"""
    calc_average_plaquette(param, Us::Array{T, 5}) where T

Calculate average plaquette for SU(2) gauge fields in 4 dimension.
"""
function calc_average_plaquette(param, Us::Array{T, 5}) where T 
    @unpack Nsite = param

    Nc = 1
    if T === SU2 
        Nc = 2
    end
    if T === SU3 
        Nc = 3
    end
    
    P = 0.0
    for ν in 1:4, μ in 1:ν-1, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        U = calc_plaquette(Us, n₁, n₂, n₃, n₄, μ, ν)
        P += 1 - real(tr(U))/Nc
    end
    P /= 6*Nsite^4
    return P
end


"""
    calc_average_plaquette(param, Us::Array{SU2, 6})

Calculate average plaquette for SU(2) gauge fields in 5 dimension.
"""
function calc_average_plaquette(param, Us::Array{SU2, 6})
    @unpack Nsite = param
    
    δ(μ, ν) = kronecker_delta(μ, ν)
    
    f(n) = check_boundary(n, Nsite)
    
    P = 0.0
    
    for ν in 1:5, μ in 1:ν-1, n₅ in 1:Nsite, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        U₁ = Us[n₁, n₂, n₃, n₄, n₅, μ]
        U₂ = Us[f(n₁+δ(1,μ)), f(n₂+δ(2,μ)), f(n₃+δ(3,μ)), f(n₄+δ(4,μ)), f(n₅+δ(5,μ)), ν]
        U₃ = conj(Us[f(n₁+δ(1,ν)), f(n₂+δ(2,ν)), f(n₃+δ(3,ν)), f(n₄+δ(4,ν)), f(n₅+δ(5,ν)), μ])
        U₄ = conj(Us[n₁, n₂, n₃, n₄, n₅, ν])
        U = U₁*U₂*U₃*U₄
        P += 1 - tr(U)/2
    end
    
    P /= 10*Nsite^5
    return P
end

