"""
    rand_U1()

Generate an element of U(1) randomly.
"""
function rand_U1()
    B = 2π*rand()
    U = exp(im*B)
end

"""
    rand_SU2()

Generate an element of SU(2) randomly.
"""
function rand_SU2()
    a₀ = 0.0 
    accept = false
    while(!accept)
        a₀ = 2*rand() - 1 # [-1,1] の一様乱数
        p  =   rand()     # [0, 1] の一様乱数
        if p ≤ sqrt(1 - a₀^2)
            accept = true
            break
        end
    end

    cosθ = 2*rand() - 1 # [-1,1] の一様乱数
    sinθ = sqrt(1 - cosθ^2)
    
    ϕ = 2π*rand() # [0, 2π] の一様乱数

    a = sqrt(1 - a₀^2)
    a₁ = a*sinθ*cos(ϕ)
    a₂ = a*sinθ*sin(ϕ)
    a₃ = a*cosθ
    
    U = SU2(a₀, a₁, a₂, a₃)
end

"""
    rand_SU3()

Generate an element of SU(3) randomly.
"""
function rand_SU3()
    U₁ = convert_SU2_to_SU3(rand_SU2(), 1)
    U₂ = convert_SU2_to_SU3(rand_SU2(), 2)
    U₃ = convert_SU2_to_SU3(rand_SU2(), 3)
    return U₃*U₂*U₁
end

"""
    initial_gaugefields!(Us::Array{ComplexF64}, param; random=false)

Generate an initial configuration of U(1) gauge fields.
"""
function initialize_gaugefields!(Us::Array{ComplexF64}, param; random=false)
    if !random
        for i in eachindex(Us)
            Us[i] = 1.0
        end
    else
        for i in eachindex(Us)
            Us[i] = rand_U1()
        end
    end
end

"""
    initial_gaugefields!(Us::Array{SU2}, param; random=false)

Generate an initial configuration of SU(2) gauge fields.
"""
function initialize_gaugefields!(Us::Array{SU2}, param; random=false)
    if !random 
        for i in eachindex(Us)
            Us[i] = SU2(1, 0, 0, 0)
        end
    else
        for i in eachindex(Us)
            Us[i] = rand_SU2()
        end
    end
end

"""
    initialize_gaugefields!(Us::Array{SU3}, param; random=false)

Generate an initial configuration of SU(3) gauge fields.
"""
function initialize_gaugefields!(Us::Array{SU3}, param; random=false)
    if !random 
        for i in eachindex(Us)
            Us[i] = SA[1 0 0;
                       0 1 0;
                       0 0 1]
        end
    else
        for i in eachindex(Us)
            Us[i] = rand_SU3()
        end
    end
end




