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
    k = Base.abs(U_staple)
        
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


"""
    generate_new_link(U_staple::SU2, U_old::SU2, β) 

Generate a new link for SU(2).
"""
function generate_new_link(U_staple::SU2, U_old::SU2, β)
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


"""
    generate_new_link(U_staple::SU3, U_old::SU3, β) 

Generate a new link for SU(3).
"""
function generate_new_link(U_staple::SU3, U_old::SU3, β)
    U = U_old

    for i in 1:3
        r = project_onto_SU2(submatrix(U*U_staple, i))
        α = generate_new_link(r, SU2(1), 2β/3)
        a = convert_SU2_to_SU3(α, i)
        U = a*U
    end

    return U
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





