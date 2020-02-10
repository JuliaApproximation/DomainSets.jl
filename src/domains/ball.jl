
# The type hierarchy is as follows:
# abstract HyperBall
# |-> abstract UnitHyperBall: radius is 1
#     |-> FixedUnitBall: dimension is part of type
#     |-> FlexibleUnitBall: dimension is specified by int field
# There are aliases for SVector{N,T} and Vector{T}:
#   EuclideanHyperBall, VectorHyperBall,
#   EuclideanUnitBall (FixedUnitBall), VectorUnitBall (FlexibleUnitBall).

"""
Supertype of balls for which elements satisfy `norm(x) < radius(ball)` (open ball)
or `norm(x) < radius(ball)` (closed ball).
"""
abstract type HyperBall{T,C} <: Domain{T} end

"A ball in a fixed N-dimensional Euclidean space."
const EuclideanHyperBall{N,T,C} = HyperBall{SVector{N,T},C}

"A ball with vector elements of variable length."
const VectorHyperBall{T,C} = HyperBall{Vector{T},C}

isclosed(::HyperBall{T,:closed}) where {T} = true
isclosed(::HyperBall{T,:open}) where {T} = false

indomain(x, ball::HyperBall) = norm(x) <= radius(ball)
approx_indomain(x, ball::HyperBall, tolerance) = norm(x) <= radius(ball)+tolerance

# A closed ball always contains at least the origin.
isempty(ball::HyperBall{T,:closed}) where {T} = false
isempty(ball::HyperBall{T,:open}) where {T} = radius(ball) == 0


"The unit ball."
abstract type UnitHyperBall{T,C} <: HyperBall{T,C} end

radius(::UnitHyperBall) = 1

"The unit ball with fixed dimension(s) specified by the element type."
struct FixedUnitBall{T,C} <: UnitHyperBall{T,C}
end

FixedUnitBall{T}() where {T} = FixedUnitBall{T,:closed}()

"The unit ball in a fixed N-dimensional space."
const EuclideanUnitBall{N,T,C} = FixedUnitBall{SVector{N,T},C}

const UnitDisk{T} = EuclideanUnitBall{2,T,:closed}
const UnitBall{T} = EuclideanUnitBall{3,T,:closed}

UnitDisk() = UnitDisk{Float64}()
UnitBall() = UnitBall{Float64}()


"The unit ball with variable dimension."
struct FlexibleUnitBall{T,C} <: UnitHyperBall{T,C}
    dimension   ::  Int
end

FlexibleUnitBall{T}(dimension::Int) where {T} =
    FlexibleUnitBall{T,:closed}(dimension)

dimension(ball::FlexibleUnitBall) = ball.dimension

"The unit ball with vector elements of a given dimension."
const VectorUnitBall{T,C} = FlexibleUnitBall{Vector{T},C}

VectorUnitBall(dimension::Int = 3) = VectorUnitBall{Float64}(dimension)

indomain(x, ball::FlexibleUnitBall) =
    (length(x) == dimension(ball)) && (norm(x) <= radius(ball))
approx_indomain(x, ball::FlexibleUnitBall, tolerance) =
    (length(x) == dimension(ball)) && (norm(x) <= radius(ball)+tolerance)

convert(::Type{Domain{T}}, ball::FixedUnitBall{S,C}) where {S,T,C} =
    FixedUnitBall{T,C}()
convert(::Type{Domain{T}}, ball::FlexibleUnitBall{S,C}) where {S,T,C} =
    FlexibleUnitBall{T,C}(ball.dimension)


show(io::IO, d::UnitHyperBall{T,:closed}) where {T} =
    print(io, "the $(dimension(d))-dimensional closed unit ball")
show(io::IO, d::UnitHyperBall{T,:open}) where {T} =
    print(io, "the $(dimension(d))-dimensional open unit ball")

# We choose the origin here
point_in_domain(::HyperBall{T}) where {T} = zero(T)



# The type hierarchy of spheres parallels that of Ball above:
# abstract HyperSphere
# |-> abstract UnitHyperSphere: radius is 1
#     |-> FixedUnitSphere: dimension is part of type
#     |-> FlexibleUnitSphere: dimension is specified by int field
# There are aliases for SVector{N,T} and Vector{T}.

"Supertype of spherical domains for which elements satisfy `norm(x) == radius(sphere)`."
abstract type HyperSphere{T} <: Domain{T} end

indomain(x, sphere::HyperSphere) = norm(x) == radius(sphere)
approx_indomain(x, sphere::HyperSphere, tolerance) =
    radius(sphere)-tolerance <= norm(x) <= radius(sphere)+tolerance

isempty(::HyperSphere) = false
isclosed(::HyperSphere) = true

"A hypersphere in a fixed N-dimensional Euclidean space."
const EuclideanHyperSphere{N,T} = HyperSphere{SVector{N,T}}

"A sphere with vector elements of variable length."
const VectorHyperSphere{T} = HyperSphere{Vector{T}}


"The unit sphere."
abstract type UnitHyperSphere{T} <: HyperSphere{T} end

radius(::UnitHyperSphere) = 1

"The unit sphere with fixed dimension(s) specified by the element type."
struct FixedUnitSphere{T} <: UnitHyperSphere{T}
end

"The unit sphere in a fixed N-dimensional Euclidean space."
const EuclideanUnitSphere{N,T} = FixedUnitSphere{SVector{N,T}}

EuclideanUnitSphere{N}() where {N} = EuclideanUnitSphere{N,Float64}()


"The unit circle in 2D."
const UnitCircle{T} = EuclideanUnitSphere{2,T}
"The unit sphere in 3D."
const UnitSphere{T} = EuclideanUnitSphere{3,T}



"The unit sphere with variable dimension."
struct FlexibleUnitSphere{T} <: UnitHyperSphere{T}
    dimension   ::  Int
end

dimension(sphere::FlexibleUnitSphere) = sphere.dimension

"The unit sphere with vector elements of a given dimension."
const VectorUnitSphere{T} = FlexibleUnitSphere{Vector{T}}

VectorUnitSphere(dimension::Int = 3) = VectorUnitSphere{Float64}(dimension)

convert(::Type{Domain{SVector{N,T}}}, d::EuclideanUnitSphere{N,S}) where {N,S,T} =
    EuclideanUnitSphere{N,T}()
convert(::Type{Domain{Vector{T}}}, d::EuclideanUnitSphere{N,S}) where {N,S,T} =
    VectorUnitSphere{T}(N)

show(io::IO, d::UnitHyperSphere) =
    dimension(d) == 2 ? print(io, "the unit circle") : print(io, "the $(dimension(d))-dimensional unit sphere")



boundary(::UnitHyperBall{N,T}) where {N,T} = UnitHyperSphere{N,T}()


################
# Applications
################

"Create a cylinder with given radius and length."
cylinder(::Type{T} = Float64) where {T} = UnitDisk{T}() × UnitInterval{T}()
cylinder(radius::T, length::T) where {T} = radius * UnitDisk{T}() × (0 .. length)

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

domain(d::UnitCircleMap{S}) where {S} = Interval{:closed,:open,S}(0, 1)

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

range(m::AngleMap{S,T}) where {S,T} = Interval{:closed,:open,T}(0, 1)

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
