# Routines having to do with computations involving domains

+(domain::Domain) = domain

function +(domain::Domain{T}, x::S) where {T,S}
    c = convert(T, x)
    (+)(domain, c)
end

*(map::Map, domain::Domain) = mapped_domain(inv(map), domain)

*(a::Number, domain::Domain{T}) where {T} = convert(Map{T}, a) * domain
*(domain::Domain, a::Number) = a*domain

/(domain::Domain{T}, a::Number) where {T} = mapped_domain(convert(Map{T}, a), domain)

+(d::Domain, x::SVector{N,T}) where {N,T} = Translation(x) * d

# Assume commutativity
+(x::AbstractVector, d::Domain) = d+x
