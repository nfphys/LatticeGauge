"""
    calc_staple(U_zero, Us, n₁, n₂, μ)

Calculate staple for link (n₁, n₂, μ) in 2 dimension.
"""
function calc_staple(U_zero, Us, n₁, n₂, μ)
    Nsite = size(Us, 1)

    δ(μ, ν) = kronecker_delta(μ, ν)
    
    f(n) = check_boundary(n, Nsite)

    U_staple = U_zero

    for ν in 1:2
        if ν === μ continue end
        
        U₁ =      Us[f(n₁+δ(1,μ)), f(n₂+δ(2,μ)), ν]
        U₂ = conj(Us[f(n₁+δ(1,ν)), f(n₂+δ(2,ν)), μ])
        U₃ = conj(Us[n₁, n₂, ν])
        U_staple += U₁*U₂*U₃
        
        U₁ = conj(Us[f(n₁+δ(1,μ)-δ(1,ν)), f(n₂+δ(2,μ)-δ(2,ν)), ν])
        U₂ = conj(Us[f(n₁-δ(1,ν)), f(n₂-δ(2,ν)), μ])
        U₃ = Us[f(n₁-δ(1,ν)), f(n₂-δ(2,ν)), ν]
        U_staple += U₁*U₂*U₃
    end
    return U_staple 
end


"""
    calc_staple(U_zero, Us, n₁, n₂, n₃, n₄, μ)

Calculate staple for link (n₁, n₂, n₃, n₄, μ) in 4 dimension.
"""
function calc_staple(U_zero, Us, n₁, n₂, n₃, n₄, μ)
    Nsite = size(Us, 1)

    δ(μ, ν) = kronecker_delta(μ, ν)
    
    f(n) = check_boundary(n, Nsite)

    U_staple = U_zero

    for ν in 1:4
        if ν === μ continue end
        
        U₁ =      Us[f(n₁+δ(1,μ)), f(n₂+δ(2,μ)), f(n₃+δ(3,μ)), f(n₄+δ(4,μ)), ν]
        U₂ = conj(Us[f(n₁+δ(1,ν)), f(n₂+δ(2,ν)), f(n₃+δ(3,ν)), f(n₄+δ(4,ν)), μ])
        U₃ = conj(Us[n₁, n₂, n₃, n₄, ν])
        U_staple += U₁*U₂*U₃
        
        U₁ = conj(Us[f(n₁+δ(1,μ)-δ(1,ν)), f(n₂+δ(2,μ)-δ(2,ν)), f(n₃+δ(3,μ)-δ(3,ν)), f(n₄+δ(4,μ)-δ(4,ν)), ν])
        U₂ = conj(Us[f(n₁-δ(1,ν)), f(n₂-δ(2,ν)), f(n₃-δ(3,ν)), f(n₄-δ(4,ν)), μ])
        U₃ = Us[f(n₁-δ(1,ν)), f(n₂-δ(2,ν)), f(n₃-δ(3,ν)), f(n₄-δ(4,ν)), ν]
        U_staple += U₁*U₂*U₃
    end
    return U_staple 
end


"""
    calc_staple(U_zero, Us, n₁, n₂, n₃, n₄, n₅, μ)

Calculate staple for link (n₁, n₂, n₃, n₄, n₅, μ) in 5 dimension.
"""
function calc_staple(U_zero, Us, n₁, n₂, n₃, n₄, n₅, μ)
    Nsite = size(Us, 1)

    δ(μ, ν) = kronecker_delta(μ, ν)
    
    f(n) = check_boundary(n, Nsite)

    U_staple = U_zero

    for ν in 1:5
        if ν === μ continue end 
        
        U₁ =      Us[f(n₁+δ(1,μ)), f(n₂+δ(2,μ)), f(n₃+δ(3,μ)), f(n₄+δ(4,μ)), f(n₅+δ(5,μ)), ν]
        U₂ = conj(Us[f(n₁+δ(1,ν)), f(n₂+δ(2,ν)), f(n₃+δ(3,ν)), f(n₄+δ(4,ν)), f(n₅+δ(5,ν)), μ])
        U₃ = conj(Us[n₁, n₂, n₃, n₄, n₅, ν])
        U_staple += U₁*U₂*U₃
        
        U₁ = conj(Us[f(n₁+δ(1,μ)-δ(1,ν)), f(n₂+δ(2,μ)-δ(2,ν)), f(n₃+δ(3,μ)-δ(3,ν)), f(n₄+δ(4,μ)-δ(4,ν)), f(n₅+δ(5,μ)-δ(5,ν)), ν])
        U₂ = conj(Us[f(n₁-δ(1,ν)), f(n₂-δ(2,ν)), f(n₃-δ(3,ν)), f(n₄-δ(4,ν)), f(n₅-δ(5,ν)), μ])
        U₃ =      Us[f(n₁-δ(1,ν)), f(n₂-δ(2,ν)), f(n₃-δ(3,ν)), f(n₄-δ(4,ν)), f(n₅-δ(5,ν)), ν]
        U_staple += U₁*U₂*U₃
    end
    return U_staple 
end