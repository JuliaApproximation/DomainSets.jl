Random.gentype(::Type{<:Domain{T}}) where T = T

Base.rand(rng::AbstractRNG, s::Random.SamplerTrivial{<:ProductDomain}) = map(i->(rand(rng, i)), factors(s[]))

Base.rand(rng::AbstractRNG, s::Random.SamplerTrivial{<:SimpleLazyDomain}) = rand(rng, superdomain(s[]))

function Base.rand(rng::AbstractRNG, s::Random.SamplerTrivial{<:Ball})
    # Technical details: http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/

    b = s[]
    
    # for low dimensions, use rejection sampling - acceptance rate is at least 52%
    if dimension(b) <= 3
        while true
            r = rand(rng, boundingbox(b))
            if r in b
                return r
            end
        end
        
    # for higher dimensions, use the "Mueller" method
    else
        u = randn_dimension(rng, eltype(b), dimension(b))
        r = rand(rng)^(1/dimension(b))
        return u.*(r/norm(u))
    end
end

randn_dimension(rng::AbstractRNG, t::Type{<:StaticVector}, d) = randn(rng, t)
randn_dimension(rng::AbstractRNG, t::Type{<:Vector}, d) = randn(rng, eltype(t), d)

# Technical notes
# ===============
#
# The methods implemented above are the only easy ones
