"""
    overrelaxation!(Us::Array{ComplexF64, 5}, param, β)

Update U(1) gauge fields by overrelaxation method in 4 dimension.
"""
function overrelaxation!(Us::Array{ComplexF64, 5}, param, β)
    return
end


"""
    overrelaxation!(Us::Array{SU2, 5}, param, β)

Update SU(2) gauge fields by overrelaxation method in 4 dimension.
"""
function overrelaxation!(Us::Array{SU2, 5}, param, β)
    @unpack Nsite = param 

    for μ in 1:4, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        
        U_zero = SU2(0,0,0,0)
        U_staple = calc_staple(U_zero, Us, n₁, n₂, n₃, n₄, μ)

        V₀ = U_staple/abs(U_staple)

        U_old = Us[n₁, n₂, n₃, n₄, μ]
        Us[n₁, n₂, n₃, n₄, μ] = V₀*U_old*conj(V₀)
    end 
end

"""
    overrelaxation!(Us::Array{SU2, 6}, param, β)

Update SU(2) gauge fields by overrelaxation method in 5 dimension.
"""
function overrelaxation!(Us::Array{SU2, 6}, param, β)
    @unpack Nsite = param 

    for μ in 1:5, n₅ in 1:Nsite, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        
        U_zero = SU2(0,0,0,0)
        U_staple = calc_staple(U_zero, Us, n₁, n₂, n₃, n₄, n₅, μ)

        V₀ = U_staple/abs(U_staple)

        U_old = Us[n₁, n₂, n₃, n₄, n₅, μ]
        Us[n₁, n₂, n₃, n₄, n₅, μ] = V₀*U_old*conj(V₀)
    end 
end


"""
    overrelaxation!(Us::Array{SU3, 5}, param, β)

Update SU(3) gauge fields by overrelaxation method in 4 dimension.
"""
function overrelaxation!(Us::Array{SU3, 5}, param, β)
    return
end