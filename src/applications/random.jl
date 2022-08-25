Random.gentype(::Type{<:Domain{T}}) where T = T

Base.rand(rng::AbstractRNG, s::Random.SamplerTrivial{<:ProductDomain}) = toexternalpoint(s[], map(i->(rand(rng, i)), factors(s[])))

Base.rand(rng::AbstractRNG, s::Random.SamplerTrivial{<:SimpleLazyDomain}) = toexternalpoint(s[], rand(rng, superdomain(s[])))

function Base.rand(rng::AbstractRNG, s::Random.SamplerTrivial{<:Ball})
    # Technical details: http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/

    b = s[]
    
    # for low-dimensional balls, use rejection sampling - acceptance rate is at least 52%
    if dimension(b) <= 3
        bb = boundingbox(b)
        while true
            r = rand(rng, bb)
            if r in b
                return r
            end
        end
        
    # for higher dimensional balls, use the "Mueller" method
    else
        u = randn_dimension(rng, eltype(b), dimension(b))
        r = radius(b)*rand(rng)^(1/dimension(b))
        return (r/norm(u))*u + center(b)
    end
end

randn_dimension(rng::AbstractRNG, t::Type{<:StaticVector}, d) = randn(rng, t)
randn_dimension(rng::AbstractRNG, t::Type{<:Vector}, d) = randn(rng, eltype(t), d)

# Implementation notes
# ====================
#
# The methods implemented above are the easy ones
# Unions and intersections could be implemented with rejection sampling, but it might be inefficient
# Sphere will require some decisions because `rand(sphere) in sphere` will usually only be approximately satisfied
# Maps may be difficult because the map could distort the distribution so that it is not uniform.
