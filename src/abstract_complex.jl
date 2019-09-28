abstract type AbstractComplex <: Number end

for f in [:exp, :sin, :cos, :tan]
    @eval Base.$f(z::AbstractComplex) = $f(Complex(z))
end
