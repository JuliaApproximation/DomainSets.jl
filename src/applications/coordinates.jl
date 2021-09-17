
"A `DomainPoint` is a point which is an element of a domain by construction."
abstract type DomainPoint{T} end

in(p::DomainPoint, d::Domain) = domain(p) == d || in(point(p), d)


## Points on a sphere

"A point on the unit sphere."
abstract type SpherePoint{T} <: DomainPoint{T} end

domain(p::SpherePoint{T}) where {T<:StaticTypes} = UnitSphere{T}()
domain(p::SpherePoint{T}) where {T<:AbstractVector} = UnitSphere{T}(length(point(p)))

"A point on the unit sphere represented by a standard Euclidean vector."
struct EuclideanSpherePoint{T} <: SpherePoint{T}
    x   ::  T
end
point(p::EuclideanSpherePoint) = p.x


"A point on the unit sphere represented in spherical coordinates."
struct SphericalCoordinate{T} <: SpherePoint{SVector{3,T}}
    θ   ::  T   # inclination or polar angle
    ϕ   ::  T   # azimuthal angle
end

point(p::SphericalCoordinate) = SVector(sin(p.θ)*cos(p.ϕ), sin(p.θ)*sin(p.ϕ), cos(p.θ))
