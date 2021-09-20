function kronecker_delta(μ, ν)
    ifelse(μ===ν, 1, 0)
end

function check_boundary(n, Nsite)
    if n < 1 
        return Nsite 
    elseif n > Nsite 
        return 1
    else
        return n
    end
end




"""
    generate_new_link(U_staple::U1, U_old::U1, β) 

Generate a new link for U(1).
"""
function generate_new_link(U_staple::U1, U_old::U1, β)
    k = abs(U_staple)
        
    accept = false
    B = 0.0
    y = 0.0
    while(!accept)
        B = 2π*rand()
        y = exp(β*k)*rand()
        if y ≤ exp(β*k*cos(B)) 
            accept = true
            break
        end
    end

    return exp(im*B)*conj(U_staple)/k
end


function SU2_update(U_staple::SU2, β)
    k = abs(U_staple) 
        
    accept = false
    temp1 = exp(β*k)
    temp2 = 1/temp1
    z = temp2
    while(!accept)
        z = (temp1 - temp2)*rand() + temp2
        y = rand()
        if y ≤ sqrt(1 - (log(z)/(β*k))^2)
            accept = true
            break
        end
    end
    
    a₀ = log(z)/(β*k)
    
    cosθ = 2*rand() - 1
    sinθ = sqrt(1 - cosθ^2)
    ϕ = 2π*rand()

    a = sqrt(1 - a₀^2)
    a₁ = a*sinθ*cos(ϕ)
    a₂ = a*sinθ*sin(ϕ)
    a₃ = a*cosθ
    
    U = SU2(a₀, a₁, a₂, a₃)

    return U*conj(U_staple)/k
end


function SU2_update_KP(U_staple::SU2, β)
    eps = 0.000000000001

    k = abs(U_staple) 
    α = β*k

    accept = false
    δ² = 0.0
    while(!accept)
        R₁ = rand() + eps 
        R₂ = rand() + eps 

        X₁ = -log(R₁)/α
        X₂ = -log(R₂)/α

        R₃ = rand()
        C = cos(2π*R₃)^2
        
        A = X₁*C 

        δ² = X₂ + A 

        R₄ = rand()
        if R₄^2 > 1 - 0.5*δ²
            continue 
        end
        accept = true
    end
    a₀ = 1 - δ²

    cosθ = 2*rand() - 1
    sinθ = sqrt(1 - cosθ^2)
    ϕ = 2π*rand()

    a = sqrt(1 - a₀^2)
    a₁ = a*sinθ*cos(ϕ)
    a₂ = a*sinθ*sin(ϕ)
    a₃ = a*cosθ
    
    U = SU2(a₀, a₁, a₂, a₃)

    return U*conj(U_staple)/k
end



"""
    generate_new_link(U_staple::SU2, U_old::SU2, β) 

Generate a new link for SU(2).
"""
function generate_new_link(U_staple::SU2, U_old::SU2, β)
    return SU2_update_KP(U_staple, β)
end


"""
    generate_new_link(U_staple::SU3, U_old::SU3, β) 

Generate a new link for SU(3).
"""
function generate_new_link(U_staple::SU3, U_old::SU3, β)
    U = U_old
    for i in 1:3
        r = project_onto_SU2(submatrix(U*U_staple, i))
        α = SU2_update_KP(r, 2β/3)
        a = convert_SU2_to_SU3(α, i)
        U = a*U
    end

    return gram_schmidt(U)
end




"""
    heatbath!(Us::Array{T, 3}, param, β) where T

Update gauge fields by heat bath method in 2 dimension.
"""
function heatbath!(Us::Array{T, 3}, param, β) where T
    @unpack Nsite, U_zero = param 
    
    for μ in 1:2, n₂ in 1:Nsite, n₁ in 1:Nsite
        U_staple = calc_staple(U_zero, Us, n₁, n₂, μ)
        U_old = Us[n₁, n₂, μ]
        Us[n₁, n₂, μ] = generate_new_link(U_staple, U_old, β)
    end
end


"""
    heatbath!(Us::Array{T, 5}, param, β) where T

Update gauge fields by heat bath method in 4 dimension.
"""
function heatbath!(Us::Array{T, 5}, param, β) where T
    @unpack Nsite, U_zero = param 
    
    for μ in 1:4, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        U_staple = calc_staple(U_zero, Us, n₁, n₂, n₃, n₄, μ)
        U_old = Us[n₁, n₂, n₃, n₄, μ]
        Us[n₁, n₂, n₃, n₄, μ] = generate_new_link(U_staple, U_old, β)
    end 
end


"""
    heatbath!(Us::Array{T, 6}, param, β) where T

Update gauge fields by heat bath method in 5 dimension.
"""
function heatbath!(Us::Array{T, 6}, param, β) where T
    @unpack Nsite, U_zero = param 
    
    for μ in 1:5, n₅ in 1:Nsite, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        U_staple = calc_staple(U_zero, Us, n₁, n₂, n₃, n₄, n₅, μ)
        U_old = Us[n₁, n₂, n₃, n₄, n₅, μ]
        Us[n₁, n₂, n₃, n₄, n₅, μ] = generate_new_link(U_staple, U_old, β)
    end 
end





#=
function SU2_update_KP2(V, β)
    eps = 0.000000000001

    ρ₀ = 0.5*real(V[1,1] + V[2,2])
    ρ₁ = 0.5*imag(V[1,2] + V[2,1])
    ρ₂ = 0.5*real(V[1,2] - V[2,1])
    ρ₃ = 0.5*imag(V[1,1] - V[2,2]) 
    ρ  = sqrt(ρ₀*ρ₀ + ρ₁*ρ₁ + ρ₂*ρ₂ + ρ₃*ρ₃)

    V₀ = inv(V/ρ)

    k  = β*ρ 

    accept = false
    δ² = 0.0
    while(!accept)
        R₀ = rand() + eps 
        R₁ = rand() + eps 

        X₀ = -log(R₀)/k 
        X₁ = -log(R₁)/k 
        
        R₂ = rand()
        C  = cos(2π*R₂)^2 

        A  = X₀*C 
        δ² = X₁ + A 
        
        R₃ = rand()
        if R₃*R₃ > 1 - 0.5*δ² 
            continue 
        end
        accept = true 
    end
    a₀ = 1 - δ² 

    rr = sqrt(1 - a₀^2)
    ϕ = 2π*rand()
    cosθ = 2*rand() - 1
    sinθ = sqrt(1 - cosθ^2)

    a₁ = rr*sinθ*cos(ϕ)
    a₂ = rr*sinθ*sin(ϕ)
    a₃ = rr*cosθ

    U_new = SA[a₀+im*a₃ im*a₁+a₂
               im*a₁-a₂ a₀-im*a₃]*V₀

    α = 0.5*(U_new[1,1] + conj(U_new[2,2]))
    β = 0.5*(U_new[2,1] - conj(U_new[1,2]))

    detU = abs2(α) + abs2(β)
    
    Unew = SA[α/detU -conj(β)/detU;
              β/detU  conj(α/detU)]
end


function project_onto_SU2!(S)
    α = 0.5*(S[1,1] + conj(S[2,2]))
    β = 0.5*(S[2,1] - conj(S[1,2]))
    SA[α -conj(β)
       β  conj(α)]
end


function make_submatrix(U, i)
    if i === 1
        return SA[U[1,1] U[1,2];
                  U[2,1] U[2,2]]
    end
    if i === 2
        return SA[U[2,2] U[2,3];
                  U[3,2] U[3,3]]
    end
    if i === 3
        return SA[U[1,1] U[1,3];
                  U[3,1] U[3,3]]
    end
end

function make_largematrix(U, i)
    a = U[1,1]
    b = U[1,2]
    c = U[2,1]
    d = U[2,2]

    if i === 1
        return SA[a b 0
                  c d 0
                  0 0 1]
    end
    if i === 2
        return SA[1 0 0 
                  0 a b
                  0 c d]
    end
    if i === 3
        return SA[a 0 b 
                  0 1 0
                  c 0 d]
    end
end
=#





