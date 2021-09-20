import Base.+
import Base.-
import Base.*
import Base./
import Base.conj
import Base.abs
import Base.abs2
import LinearAlgebra.tr
import LinearAlgebra.inv

export SU2, SU2_zero

"""
    SU2

type of an element of SU(2).
"""
struct SU2
    a₀::Float64
    a₁::Float64
    a₂::Float64 
    a₃::Float64 
end

SU2(x) = SU2(x,0,0,0)
SU2_zero() = SU2(0, 0, 0, 0)


function *(U::SU2, V::SU2)
    a₀ = U.a₀*V.a₀ - (U.a₁*V.a₁ + U.a₂*V.a₂ + U.a₃*V.a₃)
    
    a₁ = U.a₁*V.a₀ + U.a₀*V.a₁ - (U.a₂*V.a₃ - U.a₃*V.a₂)
    a₂ = U.a₂*V.a₀ + U.a₀*V.a₂ - (U.a₃*V.a₁ - U.a₁*V.a₃)
    a₃ = U.a₃*V.a₀ + U.a₀*V.a₃ - (U.a₁*V.a₂ - U.a₂*V.a₁)
    
    SU2(a₀, a₁, a₂, a₃)
end

*(U::SU2, x) = SU2(U.a₀*x, U.a₁*x, U.a₂*x, U.a₃*x)
*(x, U::SU2) = *(U, x)
/(U::SU2, x) = *(U, 1/x)

function +(U::SU2, V::SU2)
    a₀ = U.a₀ + V.a₀
    a₁ = U.a₁ + V.a₁
    a₂ = U.a₂ + V.a₂
    a₃ = U.a₃ + V.a₃
    
    SU2(a₀, a₁, a₂, a₃)
end

function -(U::SU2, V::SU2)
    a₀ = U.a₀ - V.a₀
    a₁ = U.a₁ - V.a₁
    a₂ = U.a₂ - V.a₂
    a₃ = U.a₃ - V.a₃
    
    SU2(a₀, a₁, a₂, a₃)
end

function abs2(U::SU2)
    return U.a₀*U.a₀ + U.a₁*U.a₁ + U.a₂*U.a₂ + U.a₃*U.a₃
end

abs(U::SU2) = sqrt(abs2(U))

conj(U::SU2) = SU2(U.a₀, -U.a₁, -U.a₂, -U.a₃)
inv(U::SU2) = conj(U::SU2)/abs(U)

tr(U::SU2) = 2U.a₀


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
