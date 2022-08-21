Random.gentype(::Type{<:Domain{T}}) where T = T

# Random.gentype(::Type{<:ProductDomain{T}}) where T = T
Base.rand(rng::AbstractRNG, s::Random.SamplerTrivial{<:ProductDomain}) = map(i->(rand(rng, i)), factors(s[]))

# Random.gentype(::Type{DerivedDomain{T}}) where T = Random.gentype(D)
Base.rand(rng::AbstractRNG, s::Random.SamplerTrivial{<:SimpleLazyDomain}) = rand(rng, superdomain(s[]))

# Ball
# Mapped


# Rectangle
# Wrapped
# Product Domain

# # Maybe
# Sphere
