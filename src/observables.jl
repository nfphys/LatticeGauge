"""
    calc_average_plaquette(param, Us::Array{ComplexF64})

Calculate average plaquette for U(1) gauge fields.
"""
function calc_average_plaquette(param, Us::Array{ComplexF64})
    @unpack Nsite = param
    
    δ(μ, ν) = kronecker_delta(μ, ν)
    
    f(n) = check_boundary(n, Nsite)
    
    P = 0.0
    
    for ν in 1:4, μ in 1:ν-1, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        U₁ = Us[n₁, n₂, n₃, n₄, μ]
        U₂ = Us[f(n₁+δ(1,μ)), f(n₂+δ(2,μ)), f(n₃+δ(3,μ)), f(n₄+δ(4,μ)), ν]
        U₃ = conj(Us[f(n₁+δ(1,ν)), f(n₂+δ(2,ν)), f(n₃+δ(3,ν)), f(n₄+δ(4,ν)), μ])
        U₄ = conj(Us[n₁, n₂, n₃, n₄, ν])
        U = U₁*U₂*U₃*U₄
        P += 1 - real(U)
    end
    
    P /= 6*Nsite^4
    return P
end

"""
    calc_average_plaquette(param, Us::Array{SU2})

Calculate average plaquette for SU(2) gauge fields.
"""
function calc_average_plaquette(param, Us::Array{SU2})
    @unpack Nsite = param
    
    δ(μ, ν) = kronecker_delta(μ, ν)
    
    f(n) = check_boundary(n, Nsite)
    
    P = 0.0
    
    for ν in 1:4, μ in 1:ν-1, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        U₁ = Us[n₁, n₂, n₃, n₄, μ]
        U₂ = Us[f(n₁+δ(1,μ)), f(n₂+δ(2,μ)), f(n₃+δ(3,μ)), f(n₄+δ(4,μ)), ν]
        U₃ = conj(Us[f(n₁+δ(1,ν)), f(n₂+δ(2,ν)), f(n₃+δ(3,ν)), f(n₄+δ(4,ν)), μ])
        U₄ = conj(Us[n₁, n₂, n₃, n₄, ν])
        U = U₁*U₂*U₃*U₄
        P += 1 - trace(U)/2
    end
    
    P /= 6*Nsite^4
    return P
end

function calc_wilson_loop(param, Us::Array{SU2}, I, J)
    @unpack Nsite = param
    
    δ(μ, ν) = kronecker_delta(μ, ν)
    
    f(n) = check_boundary(n, Nsite)
    
    if I < 0 || J < 0
        return NaN
    end
    
    if I === 0 || J === 0
        return 1.0
    end
    
    W = 0.0
    
    for ν in 1:4, μ in 1:ν-1, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        U = SU2(1, 0, 0, 0)
        for i in 0:I-1
            U *= Us[f(n₁+i*δ(1,μ)), f(n₂+i*δ(2,μ)), f(n₃+i*δ(3,μ)), f(n₄+i*δ(4,μ)), μ]
        end
        
        for j in 0:J-1
            U *= Us[f(n₁+I*δ(1,μ)+j*δ(1,ν)), f(n₂+I*δ(2,μ)+j*δ(2,ν)), f(n₃+I*δ(3,μ)+j*δ(3,ν)), f(n₄+I*δ(4,μ)+j*δ(4,ν)), ν]
        end
        
        for i in I-1: -1: 0
            U *= conj(Us[f(n₁+i*δ(1,μ)+J*δ(1,ν)), f(n₂+i*δ(2,μ)+J*δ(2,ν)), f(n₃+i*δ(3,μ)+J*δ(3,ν)), f(n₄+i*δ(4,μ)+J*δ(4,ν)), μ])
        end
        
        for j in J-1: -1: 0
            U *= conj(Us[f(n₁+j*δ(1,ν)), f(n₂+j*δ(2,ν)), f(n₃+j*δ(3,ν)), f(n₄+j*δ(4,ν)), ν])
        end
        W += 0.5*trace(U)
    end
    
    W /= 6*Nsite^4
end