struct ComplexPolar{R<:Real,Φ<:AbstractFloat} <: Number
    r::R
    ϕ::Φ
    function ComplexPolar(r::R, ϕ::Φ) where
        {R<:Real,Φ<:AbstractFloat}
        r >= 0 || throw(ArgumentError("r must be >= 0"))
        iszero(r) && return new{R,Φ}(r, zero(ϕ))
        return new{R,Φ}(r, mod2pi(ϕ))
    end
end

ComplexPolar(r::Real, ϕ::Real) = ComplexPolar(r, float(ϕ))
ComplexPolar{R,Φ}(r::Real, ϕ::Real) where {R,Φ} = ComplexPolar(R(r), Φ(ϕ))
ComplexPolar{R}(r::Real, ϕ::Real) where {R} = ComplexPolar(R(r), ϕ)
ComplexPolar(z::ComplexPolar) = z
(CP::Type{<:ComplexPolar})(z::Number) = CP(abs(z), angle(z))

const ComplexUnit{T} = ComplexPolar{One,T}
ComplexUnit(ϕ::Real) = ComplexPolar(One(), float(ϕ))
ComplexUnit{T}(ϕ::Real) where T<:AbstractFloat = ComplexPolar(One(), T(ϕ))

function Base.show(io::IO, z::ComplexUnit{T}) where T
    print(io, "ComplexUnit{$T}($(z.ϕ))")
end

Base.abs(z::ComplexPolar) = float(z.r)
Base.angle(z::ComplexPolar) = z.ϕ

Base.real(z::ComplexPolar) = z.r * cos(z.ϕ)
Base.imag(z::ComplexPolar) = z.r * sin(z.ϕ)

Base.abs2(z::ComplexPolar) = abs(z)^2
Base.conj(z::ComplexPolar) = ComplexPolar(z.r, -z.ϕ)

Base.Complex(z::ComplexPolar) = complex(reim(z)...)

Base.promote_rule(::Type{ComplexPolar{R,Φ}}, T::Type{<:Real}) where {R,Φ} =
    ComplexPolar{promote_type(R,T),Φ}
Base.promote_rule(::Type{ComplexPolar{R,Φ}}, ::Type{Complex{T}}) where {R,Φ,T} =
    ComplexPolar{promote_type(R,T),promote_type(Φ,T)}
Base.promote_rule(::Type{ComplexPolar{R1,Φ1}}, ::Type{ComplexPolar{R2,Φ2}}) where {R1,Φ1,R2,Φ2} =
    ComplexPolar{promote_type(R1,R2),promote_type(Φ1,Φ2)}

Base.:(==)(x::ComplexPolar, y::ComplexPolar) = x.r == y.r && x.ϕ == y.ϕ
Base.:(==)(x::ComplexPolar, y::Real) = x.r == y && iszero(x.ϕ)
Base.:(==)(x::ComplexPolar, y::Complex) = Complex(x) == y
Base.:(==)(x::Number, y::ComplexPolar) = y == x

Base.isapprox(x::ComplexPolar, y::ComplexPolar; kw...) =
    isapprox(x.r, y.r; kw...) && isapprox(x.ϕ, y.ϕ; kw...)
Base.isapprox(x::ComplexPolar, y::Real; kw...) =
    isapprox(x.r, y; kw...) && isapprox(x.ϕ, zero(y); kw...)
Base.isapprox(x::ComplexPolar, y::Complex; kw...) = isapprox(Complex(x), y; kw...)
Base.isapprox(x::Number, y::ComplexPolar; kw...) = isapprox(y, x; kw...)

Base.:+(z::ComplexPolar) = z
Base.:+(z::ComplexPolar...) = ComplexPolar(complex(sum(real, z), sum(imag, z)))
Base.:+(x::ComplexPolar, y::Number) = ComplexPolar(complex(sum(real, (x, y)), sum(imag, (x, y))))
Base.:+(x::Number, y::ComplexPolar) = y + x

Base.:-(z::ComplexPolar) = ComplexPolar(z.r, -z.ϕ)
Base.:-(z::ComplexPolar...) =
    ComplexPolar(complex(foldl((x,y) -> x-real(y), z), foldl((x,y) -> x-imag(y), z)))
Base.:-(x::ComplexPolar, y::Number) = ComplexPolar(complex(real(x)-real(y)), imag(x)-imag(y))
Base.:-(x::Number, y::ComplexPolar) = ComplexPolar(complex(real(x)-real(y)), imag(x)-imag(y))

Base.:*(x::ComplexPolar, y::ComplexPolar) = ComplexPolar(x.r * y.r, x.ϕ + y.ϕ)
Base.:*(x::ComplexPolar, y::Real) = signbit(y) ?
    ComplexPolar(x.r * -y, x.ϕ + π) :
    ComplexPolar(x.r * y, x.ϕ)
Base.:*(x::Number, y::ComplexPolar) = y * x

Base.inv(z::ComplexPolar) = ComplexPolar(inv(z.r), -z.ϕ)
Base.:/(x::ComplexPolar, y::ComplexPolar) = ComplexPolar(x.r / y.r, x.ϕ - y.ϕ)
Base.:/(x::ComplexPolar, y::Real) = signbit(y) ?
    ComplexPolar(x.r / -y, x.ϕ + π) :
    ComplexPolar(x.r / y, x.ϕ)
Base.:/(x::Real, y::ComplexPolar) = signbit(x) ?
    ComplexPolar(-x / y.r, -y.ϕ + π) :
    ComplexPolar(x / y.r, -y.ϕ)

Base.:^(x::ComplexPolar, y::Real) = ComplexPolar(x.r^y, x.ϕ*y)
Base.:^(x::ComplexPolar, y::Integer) = ComplexPolar(x.r^y, x.ϕ*y)
Base.:^(x::ComplexPolar, y::Complex) =
    ComplexPolar(x.r^real(y) * exp(-imag(y)*x.ϕ), log(x.r)*imag(y) + real(y)*x.ϕ)
Base.:^(x::Number, y::ComplexPolar) = ComplexPolar(x^Complex(y))

Base.log(C::Type{<:Complex}, z::ComplexPolar) = C(log(z.r), z.ϕ)
Base.log(z::ComplexPolar) = ComplexPolar(log(Complex, z))

for f in [:exp, :sin, :cos, :tan]
    @eval Base.$f(z::ComplexPolar) = $f(Complex(z))
end
