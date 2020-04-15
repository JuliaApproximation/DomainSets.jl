# Routines having to do with computations involving domains

+(domain::Domain) = domain

function +(domain::Domain{T}, x::S) where {T,S}
    c = convert(T, x)
    (+)(domain, c)
end

*(map::AbstractMap, domain::Domain) = forwardmap_domain(map, domain)

*(a::Number, domain::Domain{T}) where {T} = LinearMap{T}(a) * domain
*(domain::Domain, a::Number) = a*domain

/(domain::Domain{T}, a::Number) where {T} = LinearMap{T}(inv(a)) * domain


+(d::Domain, x::SVector{N,T}) where {N,T} = Translation(x) * d

# Assume commutativity
+(x::AbstractVector, d::Domain) = d+x

# Rotation around the origin
rotate(d::EuclideanDomain{2}, θ) = rotation_map(θ) * d

rotate(d::EuclideanDomain{3}, phi, theta, psi) = rotation_map(phi,theta,psi) * d
# Rotation around a fixed center.
rotate(d::EuclideanDomain{2}, θ, center::SVector{T}) where {T} = Translation(center) * (rotation_map(θ) * (Translation(-center) * d))

rotate(d::EuclideanDomain{3}, phi, theta, psi, center::SVector{T}) where {T} = Translation(center) * (rotation_map(phi,theta,psi) * (Translation(-center) * d))
