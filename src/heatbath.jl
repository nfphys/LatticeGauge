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
    heatbath!(Us::Array{SU2, 5}, param, β)

Update U(1) gauge fields by heat bath method in 4 dimension.
"""
function heatbath!(Us::Array{ComplexF64, 5}, param, β)
    @unpack Nsite = param 

    δ(μ, ν) = kronecker_delta(μ, ν)
    
    f(n) = check_boundary(n, Nsite)
    
    for μ in 1:4, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        
        U_staple = 0.0 + 0.0im
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
        
        Us[n₁, n₂, n₃, n₄, μ] = exp(im*B)*conj(U_staple)/k
    end
end


"""
    heatbath!(Us::Array{SU2, 5}, param, β)

Update SU(2) gauge fields by heat bath method in 4 dimension.
"""
function heatbath!(Us::Array{SU2, 5}, param, β)
    @unpack Nsite = param 

    δ(μ, ν) = kronecker_delta(μ, ν)
    
    f(n) = check_boundary(n, Nsite)
    
    for μ in 1:4, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        
        U_staple = SU2(0,0,0,0)
        for ν in 1:4
            if ν === μ continue end 
            
            U₁ =      Us[f(n₁+δ(1,μ)), f(n₂+δ(2,μ)), f(n₃+δ(3,μ)), f(n₄+δ(4,μ)), ν]
            U₂ = conj(Us[f(n₁+δ(1,ν)), f(n₂+δ(2,ν)), f(n₃+δ(3,ν)), f(n₄+δ(4,ν)), μ])
            U₃ = conj(Us[n₁, n₂, n₃, n₄, ν])
            U_staple += U₁*U₂*U₃
            
            U₁ = conj(Us[f(n₁+δ(1,μ)-δ(1,ν)), f(n₂+δ(2,μ)-δ(2,ν)), f(n₃+δ(3,μ)-δ(3,ν)), f(n₄+δ(4,μ)-δ(4,ν)), ν])
            U₂ = conj(Us[f(n₁-δ(1,ν)), f(n₂-δ(2,ν)), f(n₃-δ(3,ν)), f(n₄-δ(4,ν)), μ])
            U₃ =      Us[f(n₁-δ(1,ν)), f(n₂-δ(2,ν)), f(n₃-δ(3,ν)), f(n₄-δ(4,ν)), ν]
            U_staple += U₁*U₂*U₃
        end
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
        
        Us[n₁, n₂, n₃, n₄, μ] = U*conj(U_staple)/k
    end 
end



"""
    heatbath!(Us::Array{SU2, 6}, param, β)

Update SU(2) gauge fields by heat bath method in 5 dimension.
"""
function heatbath!(Us::Array{SU2, 6}, param, β)
    @unpack Nsite = param 

    δ(μ, ν) = kronecker_delta(μ, ν)
    
    f(n) = check_boundary(n, Nsite)
    
    for μ in 1:5, n₅ in 1:Nsite, n₄ in 1:Nsite, n₃ in 1:Nsite, n₂ in 1:Nsite, n₁ in 1:Nsite
        
        U_staple = SU2(0,0,0,0)
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
        
        Us[n₁, n₂, n₃, n₄, n₅, μ] = U*conj(U_staple)/k
    end 
end