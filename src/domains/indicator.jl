
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
convert(::Type{IndicatorFunction}, d::Domain{T}) where {T} =
    IndicatorFunction{T}(indicatorfunction(d))

intersectdomain1(d1::IndicatorFunction, d2) = BoxedIndicatorFunction(d1.f, d2)
intersectdomain2(d1, d2::IndicatorFunction) = BoxedIndicatorFunction(d2.f, d1)

"An indicator function with a known bounding domain."
struct BoxedIndicatorFunction{F,D,T} <: AbstractIndicatorFunction{T}
    f       ::  F
    domain  ::  D
end

BoxedIndicatorFunction(f, a, b) = BoxedIndicatorFunction(f, Rectangle(a, b))
BoxedIndicatorFunction(f::F, domain::D) where {F,T,D<:Domain{T}} =
    BoxedIndicatorFunction{F,D,T}(f, domain)

indicatorfunction(d::BoxedIndicatorFunction) = d.f

domain(d::BoxedIndicatorFunction) = d.domain

indomain(x, d::BoxedIndicatorFunction) = x ∈ domain(d) && d.f(x)


similardomain(d::BoxedIndicatorFunction, ::Type{T}) where {T} =
    BoxedIndicatorFunction(d.f, convert(Domain{T}, d.domain))

Domain(gen::Base.Generator) = generator_domain(gen)

generator_domain(gen::Base.Generator{<:Domain}) = BoxedIndicatorFunction(gen.f, gen.iter)
generator_domain(gen::Base.Generator{<:Base.Iterators.ProductIterator}) =
    productgenerator_domain(gen, gen.iter.iterators)

function productgenerator_domain(gen, domains::Tuple{Vararg{Domain,N} where N})
    domain = TupleProductDomain(gen.iter.iterators)
    BoxedIndicatorFunction(gen.f, domain)
end

boundingbox(d::BoxedIndicatorFunction) = boundingbox(d.domain)

show(io::IO, d::BoxedIndicatorFunction) =
    print(io, "indicator function bounded by: $(domain(d))")
