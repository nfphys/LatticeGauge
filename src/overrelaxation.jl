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
    @unpack Nsite = param 
    #=
    for μ in 1:4, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        
        U_zero = SU3_zero()
        U_staple = calc_staple(U_zero, Us, n₁, n₂, n₃, n₄, μ)

        v₁ = U_staple[:,1]
        v₂ = U_staple[:,2]
        v₃ = U_staple[:,3]

        v₁ /= norm(v₁)

        v₂ = v₂ - v₁*(v₁'*v₂)
        v₂ /= norm(v₂)

        v₃ = v₃ - v₁*(v₁'*v₃) - v₂*(v₂'*v₃)
        v₃ /= norm(v₃)

        V₀ = SA[v₁[1] v₂[1] v₃[1];
                v₁[2] v₂[2] v₃[2];
                v₁[3] v₂[3] v₃[3]]

        U_old = Us[n₁, n₂, n₃, n₄, μ]
        U_new = V₀*conj(U_old)*conj(V₀)

        Δ = -β/3 * real(tr((U_new - U_old)*U_staple)) 
        P = exp(-Δ)
        # Metropolis check 
        if P > 1
            Us[n₁, n₂, n₃, n₄, μ] = U_new 
        elseif rand() < P 
            Us[n₁, n₂, n₃, n₄, μ] = U_new 
        end
    end 
    =#
end