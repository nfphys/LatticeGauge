
export SU3, SU3_zero

"""
alias of SMatrix{3,3,ComplexF64, 9}.
"""
const SU3 = SMatrix{3,3,ComplexF64, 9}

"""
    SU3_identity()

Return Identity of SU(3).
"""
SU3_identity() = SA[1.0+0.0im  0.0+0.0im  0.0+0.0im;
                    0.0+0.0im  1.0+0.0im  0.0+0.0im;
                    0.0+0.0im  0.0+0.0im  1.0+0.0im]

"""
    SU3_zero()

Return zero-matrix of the same size with SU(3).
"""                
SU3_zero() = zeros(SMatrix{3,3, ComplexF64})



"""
    convert_SU2_to_SU3(U::SU2, k)

Convert an element of SU(2) to an element of SU(2) subgroups of SU(3).
"""
function convert_SU2_to_SU3(U::SU2, k)
    a = U.a₀ + im*U.a₃
    b = im*U.a₁ + U.a₂
    c = im*U.a₁ - U.a₂
    d = U.a₀ - im*U.a₃
    if k === 1
        return SA[a b 0; 
                  c d 0; 
                  0 0 1]
    end
    if k === 2
        return SA[1 0 0;
                  0 a b;
                  0 c d]
    end
    if k === 3
        return SA[d 0 c;
                  0 1 0;
                  b 0 a]
    end
end


"""
    submatrix(U::SU3, k)

Extract a 2×2 submatrix from an element of SU(3).
"""
function submatrix(U::SU3, k)
    if k === 1
        return SA[U[1,1] U[1,2];
                  U[2,1] U[2,2]]
    end
    if k === 2
        return SA[U[2,2] U[2,3];
                  U[3,2] U[3,3]]
    end
    if k === 3
        return SA[U[3,3] U[3,1];
                  U[1,3] U[1,1]]
    end
end

"""
    project_onto_SU2(U)

Project a general 2×2 matrix onto SU(2) up to normalization.
"""
function project_onto_SU2(U) 
    a = U[1,1]
    b = U[1,2]
    c = U[2,1]
    d = U[2,2]

    a₀ = real(0.5(a + d))
    a₃ = imag(0.5(a - d))

    a₁ = imag(0.5(b + c))
    a₂ = real(0.5(b - c))

    SU2(a₀, a₁, a₂, a₃)
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

