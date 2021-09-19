
export U1, U1_zero

"""
alias of ComplexF64
"""
const U1 = ComplexF64 

U1_zero() = 0.0 + 0.0im

"""
    rand_U1()

Generate an element of U(1) randomly.
"""
function rand_U1()
    B = 2π*rand()
    U = exp(im*B)
end