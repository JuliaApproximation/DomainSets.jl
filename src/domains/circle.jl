# circle.jl

"""
The unit sphere (of radius 1) in `N` dimensions.
"""
struct UnitSphere{N,T} <: EuclideanDomain{N,T}
end

const Circle{T} = UnitSphere{2,T}
const Sphere{T} = UnitBall{3,T}

indomain(x, ::UnitSphere) = norm(x) == 1

∂(::UnitBall{N,T}) where {N,T} = UnitSphere{N,T}()

circle(::Type{T} = Float64) where {T} = Circle{T}()
circle(radius::Number) = radius * circle(float(typeof(radius)))
circle(radius::Number, center::AbstractVector) = circle(radius) + center

sphere(::Type{T} = Float64) where {T} = Sphere{T}()
sphere(radius::Number) = radius * sphere(float(typeof(radius)))
sphere(radius::Number, center::AbstractVector) = sphere(radius) + center


"""
The map `[cos(2πt), sin(2πt)]` from `[0,1)` to the unit circle in `ℝ^2`.
"""
struct CircleMap{T,S} <: AbstractMap{T,S}
end

parameterization(d::Circle{T}) where {T} = CircleMap{T,subeltype(d)}()

domain(d::CircleMap{T,S}) where {T,S} = HalfOpenRightInterval{T}()

range(m::CircleMap{T,S}) where {T,S} = Circle{T}()

applymap(m::CircleMap{T,S}, t) where {T,S} = SVector(cos(2*S(pi)*t), sin(2*S(pi)*t))

function gradient(m::CircleMap{T,S}, t) where {T,S}
    a = 2*S(pi)
    SVector(-a*sin(a*t), a*cos(a*t))
end
