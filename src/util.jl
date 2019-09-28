#####
##### `One`
#####

struct One <: Real end

function One(x)
    isone(x) || throw(ArgumentError("x is not one"))
    return One()
end
One(::One) = One()

Base.promote_rule(::Type{One}, T::Type{<:Number}) = T
(T::Type{<:Number})(::One) = one(T)

Base.:*(::One, ::One) = One()
Base.:*(::One, x::Number) = x
Base.:*(x::Number, ::One) = x

Base.:/(::One, ::One) = One()
Base.:/(::One, x::Number) = inv(x)
Base.:/(x::Number, ::One) = x

Base.:^(::One, ::Number) = One()
Base.:^(::One, ::Integer) = One()

Base.log(::One) = false

Base.:(==)(::One, ::One) = true
Base.:(==)(::One, x) = isone(x)
Base.:(==)(x, ::One) = isone(x)

Base.iszero(::One) = false
Base.isone(::One) = true

Base.isapprox(::One, ::One; kw...) = true
Base.isapprox(::One, x::Number; kw...) = isapprox(x, one(x); kw...)
Base.isapprox(x::Number, ::One; kw...) = isapprox(x, one(x); kw...)
