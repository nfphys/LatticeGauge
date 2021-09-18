import Base.+
import Base.-
import Base.*
import Base./
import Base.conj

export SU2

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

trace(U::SU2) = 2U.a₀

