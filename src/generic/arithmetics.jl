# arithmetics.jl
# Routines having to do with computations involving domains

+(domain::Domain) = domain

function +(domain::Domain{T}, x::S) where {T,S}
    c = convert_space(spacetype(T), x)
    c isa S && error("Cannot add $x to $domain")
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
rotate(d::Domain2d, θ) = rotation_map(θ) * d

rotate(d::Domain3d, phi, theta, psi) = rotation_map(phi,theta,psi) * d
# Rotation around a fixed center.
rotate(d::Domain2d, θ, center::SVector{T}) where {T} = Translation(center) * (rotation_map(θ) * (Translation(-center) * d))

rotate(d::Domain3d, phi, theta, psi, center::SVector{T}) where {T} = Translation(center) * (rotation_map(phi,theta,psi) * (Translation(-center) * d))
