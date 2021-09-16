
"The set of all natural numbers."
struct NaturalNumbers <: Domain{Int}
end

in(x, d::NaturalNumbers) = isinteger(x) && (x >= 0)
in(x::AbstractArray, d::NaturalNumbers) = false
==(::NaturalNumbers, ::NaturalNumbers) = true

approx_in(x::Real, d::NaturalNumbers, tol) =
    (abs(x-round(x)) < tol) && (round(Int, x) ∈ d)


"The set of all integers."
struct Integers <: Domain{Int}
end

in(x, d::Integers) = isinteger(x)
in(x::AbstractArray, d::Integers) = false
==(::Integers, ::Integers) = true

approx_in(x::Real, d::Integers, tol) =
    (abs(x-round(x)) < tol) && (round(Int, x) ∈ d)


"The set of all real numbers."
struct RealNumbers <: Domain{Float64}
end

in(x::Number, d::RealNumbers) = isreal(x)
# isreal also allows real arrays, we want to disallow that here:
in(x::AbstractArray, d::RealNumbers) = false

approx_in(x::Complex, d::RealNumbers, tol) = (imag(x) < tol) && (real(x) ∈ d)
==(::RealNumbers, ::RealNumbers) = true

"The set of all rationals."
struct Rationals <: Domain{Rational{Int}}
end

in(x, d::Rationals) = x ∈ Integers()
in(x::Rational, d::Rationals) = true
==(::Rationals, ::Rationals) = true

"The set of all complex numbers whose real and imaginary parts are real numbers."
struct ComplexNumbers <: Domain{Complex{Float64}}
end

in(x::Complex{T}, d::ComplexNumbers) where {T} =
    isreal(real(x)) && isreal(imag(x))
in(x::Complex{T}, d::ComplexNumbers) where {T<:Real} = true
in(x, d::ComplexNumbers) = x ∈ RealNumbers()
==(::ComplexNumbers, ::ComplexNumbers) = true


"The set of natural numbers."
const ℕ = NaturalNumbers()
"The set of integers."
const ℤ = Integers()
"The set of rational numbers."
const ℚ = Rationals()
"The set of real numbers."
const ℝ = RealNumbers()
"The set of complex numbers."
const ℂ = ComplexNumbers()

"The space ℝ^1."
const ℝ1 = VcatDomain(ℝ)
"The space ℝ^2."
const ℝ2 = VcatDomain(ℝ, ℝ)
"The space ℝ^3."
const ℝ3 = VcatDomain(ℝ, ℝ, ℝ)
"The space ℝ^4."
const ℝ4 = VcatDomain(ℝ, ℝ, ℝ, ℝ)
