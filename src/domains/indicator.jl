
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
indicatorfunction(d) = x -> x ∈ d

indomain(x, d::AbstractIndicatorFunction) = _indomain(x, d, indicatorfunction(d))
_indomain(x, d::AbstractIndicatorFunction, f) = f(x)

show(io::IO, d::AbstractIndicatorFunction) =
    print(io, "indicator domain defined by function f = $(indicatorfunction(d))")


"An `IndicatorFunction` is a domain that implements `f(x) = x ∈ D` by storing `f`."
struct IndicatorFunction{T,F} <: AbstractIndicatorFunction{T}
    f       ::  F
end

IndicatorFunction(f) = IndicatorFunction{Float64}(f)
IndicatorFunction{T}(f::F) where {T,F} = IndicatorFunction{T,F}(f)

indicatorfunction(d::IndicatorFunction) = d.f

similardomain(d::IndicatorFunction, ::Type{T}) where {T} = IndicatorFunction{T}(d.f)

convert(::Type{IndicatorFunction}, d::AbstractIndicatorFunction) = d
convert(::Type{IndicatorFunction}, d) =
    IndicatorFunction{domaineltype(d)}(indicatorfunction(checkdomain(d)))

isequaldomain(d1::IndicatorFunction, d2::IndicatorFunction) =
    indicatorfunction(d1)==indicatorfunction(d2)
domainhash(d1::IndicatorFunction, h::UInt) = hashrec("IndicatorFunction", indicatorfunction(d1), h)

intersectdomain1(d1::IndicatorFunction, d2) = BoundedIndicatorFunction(d1.f, d2)
intersectdomain2(d1, d2::IndicatorFunction) = BoundedIndicatorFunction(d2.f, d1)

"An indicator function with a known bounding domain."
struct BoundedIndicatorFunction{T,F,D} <: AbstractIndicatorFunction{T}
    f       ::  F
    domain  ::  D
end

BoundedIndicatorFunction(f, domain) =
    BoundedIndicatorFunction{domaineltype(domain)}(f, domain)
BoundedIndicatorFunction{T}(f, domain::Domain{T}) where {T} =
    BoundedIndicatorFunction{T,typeof(f),typeof(domain)}(f, domain)
BoundedIndicatorFunction{T}(f, domain) where {T} =
    _BoundedIndicatorFunction(T, f, checkdomain(domain))
_BoundedIndicatorFunction(::Type{T}, f, domain) where {T} =
    BoundedIndicatorFunction{T,typeof(f),typeof(domain)}(f, domain)

indicatorfunction(d::BoundedIndicatorFunction) = d.f

boundingdomain(d::BoundedIndicatorFunction) = d.domain

indomain(x, d::BoundedIndicatorFunction) = in(x, boundingdomain(d)) && d.f(x)

isequaldomain(d1::BoundedIndicatorFunction, d2::BoundedIndicatorFunction) =
    indicatorfunction(d1)==indicatorfunction(d2) && boundingdomain(d1)==boundingdomain(d2)
domainhash(d::BoundedIndicatorFunction, h::UInt) =
    hashrec("BoundedIndicatorFunction", indicatorfunction(d), boundingdomain(d), h)

similardomain(d::BoundedIndicatorFunction, ::Type{T}) where {T} =
    BoundedIndicatorFunction(d.f, convert_eltype(T, d.domain))

boundingbox(d::BoundedIndicatorFunction) = boundingbox(boundingdomain(d))

function show(io::IO, d::BoundedIndicatorFunction)
    print(io, "indicator function bounded by: ")
    show(io, boundingdomain(d))
end
