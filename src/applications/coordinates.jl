
"""
A `DomainPoint` is a point which is an element of a domain by construction.

A domain point is just a point, not a domain. This is different from a `Point`
type.

The type is a wrapper: retrieve the underlying point using `point(p)`.
"""
abstract type DomainPoint{T} end

## Points on a sphere

"A point on the unit sphere."
abstract type SpherePoint{T} <: DomainPoint{T} end

domain(p::SpherePoint{T}) where {T<:StaticTypes} = UnitSphere{T}()
domain(p::SpherePoint{T}) where {T<:AbstractVector} = UnitSphere{T}(length(point(p)))

"A point on the unit sphere represented by a standard Euclidean vector."
struct EuclideanSpherePoint{T} <: SpherePoint{T}
    x   ::  T

    function EuclideanSpherePoint{T}(x::T) where T
        @assert norm(x) ≈ 1
        new(x)
    end
end
EuclideanSpherePoint(x::T) where T = EuclideanSpherePoint{T}(x)

point(p::EuclideanSpherePoint) = p.x

in(p::EuclideanSpherePoint{T}, d::UnitSphere{T}) where {T<:StaticTypes} = true
in(p::EuclideanSpherePoint{T}, d::UnitSphere{T}) where {T<:AbstractVector} =
    length(point(p)) == dimension(d)

"""
A point on the unit sphere represented in spherical coordinates.
"""
struct SphericalCoordinate{T} <: SpherePoint{SVector{3,T}}
    θ   ::  T   # inclination or polar angle
    ϕ   ::  T   # azimuthal angle
end
SphericalCoordinate(θ, ϕ) = SphericalCoordinate(promote(θ, ϕ)...)
SphericalCoordinate(θ::T, ϕ::T) where {T <: Integer} =
    SphericalCoordinate(float(θ), float(ϕ))

point(p::SphericalCoordinate) = SVector(sin(p.θ)*cos(p.ϕ), sin(p.θ)*sin(p.ϕ), cos(p.θ))

in(p::SphericalCoordinate{T}, d::EuclideanUnitSphere{3,T}) where T = true
