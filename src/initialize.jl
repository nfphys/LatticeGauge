

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




