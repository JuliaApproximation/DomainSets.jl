# circle.jl

abstract type AbstractHyperSphere{N,T} <: EuclideanDomain{N,T} end



"""
The unit sphere (of radius 1) in `N` dimensions.
"""
struct UnitHyperSphere{N,T} <: AbstractHyperSphere{N,T} end

UnitHyperSphere{N}() where N = UnitHyperSphere{N,Float64}()

const UnitCircle{T} = UnitHyperSphere{2,T}
const UnitSphere{T} = UnitHyperSphere{3,T}

convert(::Type{Domain{SVector{N,T}}}, d::UnitHyperSphere{N}) where {N,T} =
    UnitHyperSphere{N,T}()

indomain(x, ::UnitHyperSphere) = norm(x) == 1

approx_indomain(x, ::UnitHyperSphere, tolerance) = 1-tolerance <= norm(x) <= 1+tolerance

boundary(::UnitHyperBall{N,T}) where {N,T} = UnitHyperSphere{N,T}()

isempty(::UnitHyperSphere) = false

"Create an ellipse curve with semi-axes lengths `a` and `b` respectively."
ellipse(a::Number, b::Number) = ellipse(promote(a,b)...)
ellipse(a::T, b::T) where {T <: Number} = scaling_map(a, b) * UnitCircle{T}()

"Create an ellipse-shaped domain with semi-axes lengths `a` and `b` respectively."
ellipse_shape(a::Number, b::Number) = ellipse_shape(promote(a,b)...)
ellipse_shape(a::T, b::T) where {T <: Number} = scaling_map(a, b) * UnitDisk{T}()


"""
The map `[cos(2πt), sin(2πt)]` from `[0,1)` to the unit circle in `ℝ^2`.
"""
struct UnitCircleMap{S,T} <: AbstractMap{S,T} end

parameterization(d::UnitCircle) = UnitCircleMap{subeltype(d),eltype(d)}()

domain(d::UnitCircleMap{S}) where S = HalfOpenRightInterval{S}(0, 1)

image(m::UnitCircleMap{S}) where S = UnitCircle{S}()

applymap(m::UnitCircleMap{S}, t) where S = SVector(cos(2*S(pi)*t), sin(2*S(pi)*t))

function gradient(m::UnitCircleMap{S}, t) where S
    a = 2*S(pi)
    SVector(-a*sin(a*t), a*cos(a*t))
end


"""
`AngleMap` is a left inverse of `UnitCircleMap`. A 2D vector `x` is projected onto
the intersection point with the unit circle of the line connecting `x` to the
origin. The angle of this point, scaled to the interval `[0,1)`, is the result.
"""
struct AngleMap{S,T} <: AbstractMap{S,T}
end

domain(d::AngleMap{S,T}) where {S,T} = FullSpace{S}()

range(m::AngleMap{S,T}) where {S,T} = HalfOpenRightInterval{T}(0, 1)

function applymap(m::AngleMap, x)
    twopi = 2*convert(codomaintype(m), pi)
    θ = atan(x[2],x[1])
    if θ < 0
        # atan2 returns an angle in (-π,π], convert to [0,2π) using periodicity.
        θ += twopi
    end
    # And divide by 2π to scale to [0,1)
    θ / twopi
end

left_inverse(m::UnitCircleMap{S,T}) where {S,T} = AngleMap{T,S}()

right_inverse(m::AngleMap{S,T}) where {S,T} = UnitCircleMap{T,S}()
