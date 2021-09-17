
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
SU3_zero() = SA[0.0+0.0im  0.0+0.0im  0.0+0.0im;
                0.0+0.0im  0.0+0.0im  0.0+0.0im;
                0.0+0.0im  0.0+0.0im  0.0+0.0im]


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
    return SU3_identity()
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
    return SA[1.0+0.0im 0.0+0.0im;
              0.0+0.0im 1.0+0.0im]
end
