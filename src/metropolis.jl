"""
    metropolis!(Us::Array{ComplexF64, 5}, param, β)

Update U(1) gauge fields by heat bath method in 4 dimension.
"""
function metropolis!(Us::Array{ComplexF64, 5}, param, β)
    @unpack Nsite, Nhit = param 

    for μ in 1:4, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        
        U_zero = 0.0 + 0.0im
        U_staple = calc_staple(U_zero, Us, n₁, n₂, n₃, n₄, μ)

        for i in 1:Nhit 
            U_old = Us[n₁, n₂, n₃, n₄, μ]
            U_new = rand_U1()*U_old 

            Δ = -β * real((U_new - U_old)*U_staple)
            P = exp(-Δ)
            # Metropolis check 
            if P > 1
                Us[n₁, n₂, n₃, n₄, μ] = U_new 
            elseif rand() < P 
                Us[n₁, n₂, n₃, n₄, μ] = U_new 
            end
        end
    end 
end


"""
    metropolis!(Us::Array{SU2, 5}, param, β)

Update SU(2) gauge fields by heat bath method in 4 dimension.
"""
function metropolis!(Us::Array{SU2, 5}, param, β)
    @unpack Nsite, Nhit = param 

    for μ in 1:4, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        
        U_zero = SU2(0,0,0,0)
        U_staple = calc_staple(U_zero, Us, n₁, n₂, n₃, n₄, μ)

        for i in 1:Nhit 
            U_old = Us[n₁, n₂, n₃, n₄, μ]
            U_new = rand_SU2()*U_old 

            Δ = -β/2 * trace((U_new - U_old)*U_staple)
            P = exp(-Δ)
            # Metropolis check 
            if P > 1
                Us[n₁, n₂, n₃, n₄, μ] = U_new 
            elseif rand() < P 
                Us[n₁, n₂, n₃, n₄, μ] = U_new 
            end
        end
    end 
end

"""
    metropolis!(Us::Array{SU2, 6}, param, β)

Update SU(2) gauge fields by heat bath method in 5 dimension.
"""
function metropolis!(Us::Array{SU2, 6}, param, β)
    @unpack Nsite, Nhit = param 

    for μ in 1:5, n₅ in 1:Nsite, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        
        U_zero = SU2(0,0,0,0)
        U_staple = calc_staple(U_zero, Us, n₁, n₂, n₃, n₄, n₅, μ)

        for i in 1:Nhit 
            U_old = Us[n₁, n₂, n₃, n₄, n₅, μ]
            U_new = rand_SU2()*U_old 

            Δ = -β/2 * trace((U_new - U_old)*U_staple)
            P = exp(-Δ)
            # Metropolis check 
            if P > 1
                Us[n₁, n₂, n₃, n₄, n₅, μ] = U_new 
            elseif rand() < P 
                Us[n₁, n₂, n₃, n₄, n₅, μ] = U_new 
            end
        end
    end 
end


"""
    metropolis!(Us::Array{SU3, 5}, param, β)

Update SU(3) gauge fields by heat bath method in 4 dimension.
"""
function metropolis!(Us::Array{SU3, 5}, param, β)
    @unpack Nsite, Nhit = param 

    for μ in 1:4, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        
        U_zero = SU3_zero()
        U_staple = calc_staple(U_zero, Us, n₁, n₂, n₃, n₄, μ)

        for i in 1:Nhit 
            U_old = Us[n₁, n₂, n₃, n₄, μ]
            U_new = rand_SU3()*U_old 

            Δ = -β/3 * real(tr((U_new-U_old)*U_staple)) 
            P = exp(-Δ)
            # Metropolis check 
            if P > 1
                Us[n₁, n₂, n₃, n₄, μ] = U_new 
            elseif rand() < P 
                Us[n₁, n₂, n₃, n₄, μ] = U_new 
            end
        end
    end 
end