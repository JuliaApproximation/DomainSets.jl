# circle.jl

"""
The unit sphere (of radius 1) in `N` dimensions.
"""
struct UnitSphere{N,T} <: EuclideanDomain{N,T}
end

const Circle{T} = UnitSphere{2,T}
const Sphere{T} = UnitSphere{3,T}

indomain(x, ::UnitSphere) = norm(x) == 1

approx_indomain(x, ::UnitSphere, tolerance) = 1-tolerance <= norm(x) <= 1+tolerance

∂(::UnitBall{N,T}) where {N,T} = UnitSphere{N,T}()

circle(::Type{T} = Float64) where {T} = Circle{T}()
circle(radius::Number) = radius * circle(float(typeof(radius)))
circle(radius::Number, center::AbstractVector) = circle(radius) + center

sphere(::Type{T} = Float64) where {T} = Sphere{T}()
sphere(radius::Number) = radius * sphere(float(typeof(radius)))
sphere(radius::Number, center::AbstractVector) = sphere(radius) + center

"Create an ellipse curve with semi-axes lengths `a` and `b` respectively."
ellipse(a::Number, b::Number) = ellipse(promote(a,b)...)
ellipse(a::T, b::T) where {T <: Number} = scaling_map(a, b) * Circle{T}()

"Create an ellipse-shaped domain with semi-axes lengths `a` and `b` respectively."
ellipse_shape(a::Number, b::Number) = ellipse_shape(promote(a,b)...)
ellipse_shape(a::T, b::T) where {T <: Number} = scaling_map(a, b) * Disk{T}()


"""
The map `[cos(2πt), sin(2πt)]` from `[0,1)` to the unit circle in `ℝ^2`.
"""
struct CircleMap{S,T} <: AbstractMap{S,T}
end

parameterization(d::Circle) = CircleMap{subeltype(d),eltype(d)}()

domain(d::CircleMap{S}) where S = HalfOpenRightInterval{S}(0, 1)

image(m::CircleMap{S}) where S = Circle{S}()

applymap(m::CircleMap{S}, t) where S = SVector(cos(2*S(pi)*t), sin(2*S(pi)*t))

function gradient(m::CircleMap{S}, t) where S
    a = 2*S(pi)
    SVector(-a*sin(a*t), a*cos(a*t))
end


"""
`AngleMap` is a left inverse of `CircleMap`. A 2D vector `x` is projected onto
the intersection point with the unit circle of the line connecting `x` to the
origin. The angle of this point, scaled to the interval `[0,1)`, is the result.
"""
struct AngleMap{S,T} <: AbstractMap{S,T}
end

domain(d::AngleMap{S,T}) where {S,T} = FullSpace{S}()

range(m::AngleMap{S,T}) where {S,T} = HalfOpenRightInterval{T}(0, 1)

function applymap(m::AngleMap, x)
    twopi = 2*convert(codomaintype(m), pi)
    θ = (VERSION < v"0.7-") ? atan2(x[2],x[1]) : atan(x[2],x[1])
    if θ < 0
        # atan2 returns an angle in (-π,π], convert to [0,2π) using periodicity.
        θ += twopi
    end
    # And divide by 2π to scale to [0,1)
    θ / twopi
end

left_inverse(m::CircleMap{S,T}) where {S,T} = AngleMap{T,S}()

right_inverse(m::AngleMap{S,T}) where {S,T} = CircleMap{T,S}()
