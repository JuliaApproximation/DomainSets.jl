
"""
Supertype of domains that are defined by an indicator function.

An indicator function is a function `f : S -> [0,1]` that indicates membership
of `x` to a domain `D` with `D ⊂ S`. The indicator function corresponds exactly
to the `in` function of a domain: `f(x) = x ∈ D`.

Concrete subtypes of `AbstractIndicatorFunction` store a representation of this
indicator function and implement `in` using that representation, rather than
implementing `in` directly.
"""
abstract type AbstractIndicatorFunction{T} <: Domain{T} end

"The indicator function of a domain is the function `f(x) = x ∈ D`."
indicatorfunction(d::Domain) = x -> x ∈ d

indomain(x, d::AbstractIndicatorFunction) = _indomain(x, d, indicatorfunction(d))
_indomain(x, d::AbstractIndicatorFunction, f) = f(x)

convert(::Type{Domain{T}}, d::AbstractIndicatorFunction{T}) where {T} = d
convert(::Type{Domain{T}}, d::AbstractIndicatorFunction{S}) where {S,T} =
    similar_indicatorfunction(d, T)

show(io::IO, d::AbstractIndicatorFunction) =
    print(io, "indicator domain defined by function f = $(indicatorfunction(d))")


"An `IndicatorFunction` is a domain that implements `f(x) = x ∈ D` by storing `f`."
struct IndicatorFunction{T,F} <: AbstractIndicatorFunction{T}
    f       ::  F
end

IndicatorFunction(f) = IndicatorFunction{Float64}(f)
IndicatorFunction{T}(f::F) where {T,F} = IndicatorFunction{T,F}(f)

indicatorfunction(d::IndicatorFunction) = d.f

similar_indicatorfunction(d::IndicatorFunction, ::Type{T}) where {T} =
    IndicatorFunction{T}(d.f)
