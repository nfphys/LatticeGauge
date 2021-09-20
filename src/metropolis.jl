"""
    generate_SU2_near_identity(k)

Generate an element of SU(2) near identity with probability 
P(h) ∝ exp(k Tr(h))
"""
function generate_SU2_near_identity(k)
    accept = false
    temp1 = exp(2k)
    temp2 = 1/temp1
    z = temp2
    while(!accept)
        z = (temp1 - temp2)*rand() + temp2
        y = rand()
        if y ≤ sqrt(1 - (log(z)/(2k))^2)
            accept = true
            break
        end
    end
    
    a₀ = log(z)/(2k)
    
    cosθ = 2*rand() - 1
    sinθ = sqrt(1 - cosθ^2)
    ϕ = 2π*rand()

    a = sqrt(1 - a₀^2)
    a₁ = a*sinθ*cos(ϕ)
    a₂ = a*sinθ*sin(ϕ)
    a₃ = a*cosθ
    
    h = SU2(a₀, a₁, a₂, a₃)
end


"""
    metropolis!(Us::Array{ComplexF64, 3}, param, β)

Update U(1) gauge fields by heat bath method in 2 dimension.
"""
function metropolis!(Us::Array{ComplexF64, 3}, param, β)
    @unpack Nsite, Nhit = param 

    for μ in 1:2, n₂ in 1:Nsite, n₁ in 1:Nsite
        
        U_zero = 0.0 + 0.0im
        U_staple = calc_staple(U_zero, Us, n₁, n₂, μ)

        for i in 1:Nhit 
            U_old = Us[n₁, n₂, μ]
            U_new = rand_U1()*U_old 

            Δ = -β * real((U_new - U_old)*U_staple)
            P = exp(-Δ)
            # Metropolis check 
            if P > 1
                Us[n₁, n₂, μ] = U_new 
            elseif rand() < P 
                Us[n₁, n₂, μ] = U_new 
            end
        end
    end 
end


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
    metropolis!(Us::Array{SU2, 3}, param, β)

Update SU(2) gauge fields by heat bath method in 2 dimension.
"""
function metropolis!(Us::Array{SU2, 3}, param, β)
    @unpack Nsite, Nhit = param 

    for μ in 1:2, n₂ in 1:Nsite, n₁ in 1:Nsite
        
        U_zero = SU2(0,0,0,0)
        U_staple = calc_staple(U_zero, Us, n₁, n₂, μ)

        for i in 1:Nhit 
            U_old = Us[n₁, n₂, μ] 

            h = generate_SU2_near_identity(2β)
            U_new = h*U_old 

            Δ = -β/2 * tr((U_new - U_old)*U_staple)
            P = exp(-Δ)
            # Metropolis check 
            if P > 1
                Us[n₁, n₂, μ] = U_new 
            elseif rand() < P 
                Us[n₁, n₂, μ] = U_new 
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

            h = generate_SU2_near_identity(2β)
            U_new = h*U_old 

            Δ = -β/2 * tr((U_new - U_old)*U_staple)
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

            Δ = -β/2 * tr((U_new - U_old)*U_staple)
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
    metropolis!(Us::Array{SU3, 3}, param, β)

Update SU(3) gauge fields by heat bath method in 2 dimension.
"""
function metropolis!(Us::Array{SU3, 3}, param, β)
    @unpack Nsite, Nhit = param 

    for μ in 1:2, n₂ in 1:Nsite, n₁ in 1:Nsite
        
        U_zero = SU3_zero()
        U_staple = calc_staple(U_zero, Us, n₁, n₂, μ)

        for i in 1:Nhit 
            for k in 1:3
                U_old = Us[n₁, n₂, μ]

                h = generate_SU2_near_identity(2β)
                U_new = convert_SU2_to_SU3(h, k)*U_old

                Δ = -β/3 * real(tr((U_new - U_old)*U_staple)) 
                P = exp(-Δ)
                # Metropolis check 
                if P > 1
                    Us[n₁, n₂, μ] = U_new 
                elseif rand() < P 
                    Us[n₁, n₂, μ] = U_new 
                end
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
            for k in 1:3
                U_old = Us[n₁, n₂, n₃, n₄, μ]

                h = generate_SU2_near_identity(2β)
                U_new = convert_SU2_to_SU3(h, k)*U_old

                # Metropolis check 
                Δ = -β/3 * real(tr((U_new - U_old)*U_staple)) 
                P = exp(-Δ)
                if P > 1
                    Us[n₁, n₂, n₃, n₄, μ] = U_new 
                elseif rand() < P 
                    Us[n₁, n₂, n₃, n₄, μ] = U_new 
                end
            end
        end
    end 
end