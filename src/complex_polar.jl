struct PolarComplex{R<:Real,Φ<:AbstractFloat} <: AbstractComplex
    r::R
    ϕ::Φ
    function PolarComplex(r::R, ϕ::Φ) where
        {R<:Real,Φ<:AbstractFloat}
        r >= 0 || throw(ArgumentError("r must be >= 0"))
        iszero(r) && return new{R,Φ}(r, zero(ϕ))
        return new{R,Φ}(r, mod2pi(ϕ))
    end
end

PolarComplex(r::Real, ϕ::Real) = PolarComplex(r, float(ϕ))
PolarComplex{R,Φ}(r::Real, ϕ::Real) where {R,Φ} = PolarComplex(R(r), Φ(ϕ))
PolarComplex{R}(r::Real, ϕ::Real) where {R} = PolarComplex(R(r), ϕ)
PolarComplex(z::PolarComplex) = z
(CP::Type{<:PolarComplex})(z::Number) = CP(abs(z), angle(z))

const UnitComplex{T} = PolarComplex{One,T}
UnitComplex(ϕ::Real) = PolarComplex(One(), float(ϕ))
UnitComplex{T}(ϕ::Real) where T<:AbstractFloat = PolarComplex(One(), T(ϕ))

function Base.show(io::IO, z::UnitComplex{T}) where T
    print(io, "UnitComplex{$T}($(z.ϕ))")
end

Base.abs(z::PolarComplex) = float(z.r)
Base.angle(z::PolarComplex) = z.ϕ

Base.real(z::PolarComplex) = z.r * cos(z.ϕ)
Base.imag(z::PolarComplex) = z.r * sin(z.ϕ)

Base.abs2(z::PolarComplex) = abs(z)^2
Base.conj(z::PolarComplex) = PolarComplex(z.r, -z.ϕ)

Base.Complex(z::PolarComplex) = complex(reim(z)...)

Base.promote_rule(::Type{PolarComplex{R,Φ}}, T::Type{<:Real}) where {R,Φ} =
    PolarComplex{promote_type(R,T),Φ}
Base.promote_rule(::Type{PolarComplex{R,Φ}}, ::Type{Complex{T}}) where {R,Φ,T} =
    PolarComplex{promote_type(R,T),promote_type(Φ,T)}
Base.promote_rule(::Type{PolarComplex{R1,Φ1}}, ::Type{PolarComplex{R2,Φ2}}) where {R1,Φ1,R2,Φ2} =
    PolarComplex{promote_type(R1,R2),promote_type(Φ1,Φ2)}

Base.:(==)(x::PolarComplex, y::PolarComplex) = x.r == y.r && x.ϕ == y.ϕ
Base.:(==)(x::PolarComplex, y::Real) = x.r == y && iszero(x.ϕ)
Base.:(==)(x::PolarComplex, y::Complex) = Complex(x) == y
Base.:(==)(x::Number, y::PolarComplex) = y == x

Base.isapprox(x::PolarComplex, y::PolarComplex; kw...) =
    isapprox(x.r, y.r; kw...) && isapprox(x.ϕ, y.ϕ; kw...)
Base.isapprox(x::PolarComplex, y::Real; kw...) =
    isapprox(x.r, y; kw...) && isapprox(x.ϕ, zero(y); kw...)
Base.isapprox(x::PolarComplex, y::Complex; kw...) = isapprox(Complex(x), y; kw...)
Base.isapprox(x::Number, y::PolarComplex; kw...) = isapprox(y, x; kw...)

Base.:+(z::PolarComplex) = z
Base.:+(z::PolarComplex...) = PolarComplex(complex(sum(real, z), sum(imag, z)))
Base.:+(x::PolarComplex, y::Number) = PolarComplex(complex(sum(real, (x, y)), sum(imag, (x, y))))
Base.:+(x::Number, y::PolarComplex) = y + x

Base.:-(z::PolarComplex) = PolarComplex(z.r, -z.ϕ)
Base.:-(z::PolarComplex...) =
    PolarComplex(complex(foldl((x,y) -> x-real(y), z), foldl((x,y) -> x-imag(y), z)))
Base.:-(x::PolarComplex, y::Number) = PolarComplex(complex(real(x)-real(y)), imag(x)-imag(y))
Base.:-(x::Number, y::PolarComplex) = PolarComplex(complex(real(x)-real(y)), imag(x)-imag(y))

Base.:*(x::PolarComplex, y::PolarComplex) = PolarComplex(x.r * y.r, x.ϕ + y.ϕ)
Base.:*(x::PolarComplex, y::Real) = signbit(y) ?
    PolarComplex(x.r * -y, x.ϕ + π) :
    PolarComplex(x.r * y, x.ϕ)
Base.:*(x::Number, y::PolarComplex) = y * x

Base.inv(z::PolarComplex) = PolarComplex(inv(z.r), -z.ϕ)
Base.:/(x::PolarComplex, y::PolarComplex) = PolarComplex(x.r / y.r, x.ϕ - y.ϕ)
Base.:/(x::PolarComplex, y::Real) = signbit(y) ?
    PolarComplex(x.r / -y, x.ϕ + π) :
    PolarComplex(x.r / y, x.ϕ)
Base.:/(x::Real, y::PolarComplex) = signbit(x) ?
    PolarComplex(-x / y.r, -y.ϕ + π) :
    PolarComplex(x / y.r, -y.ϕ)

Base.:^(x::PolarComplex, y::Real) = PolarComplex(x.r^y, x.ϕ*y)
Base.:^(x::PolarComplex, y::Integer) = PolarComplex(x.r^y, x.ϕ*y)
Base.:^(x::PolarComplex, y::Complex) =
    PolarComplex(x.r^real(y) * exp(-imag(y)*x.ϕ), log(x.r)*imag(y) + real(y)*x.ϕ)
Base.:^(x::Number, y::PolarComplex) = PolarComplex(x^Complex(y))

Base.log(C::Type{<:Complex}, z::PolarComplex) = C(log(z.r), z.ϕ)
Base.log(z::PolarComplex) = PolarComplex(log(Complex, z))
