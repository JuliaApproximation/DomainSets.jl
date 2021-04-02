
# The type hierarchy is as follows:
# abstract HyperBall
# |-> abstract UnitHyperBall: radius is 1
#     |-> StaticUnitBall: dimension is part of type
#     |-> DynamicUnitBall: dimension is specified by int field
# There are aliases for SVector{N,T} and Vector{T}:
#   EuclideanHyperBall, VectorHyperBall,
#   EuclideanUnitBall (StaticUnitBall), VectorUnitBall (DynamicUnitBall).

"""
Supertype of balls for which elements satisfy `norm(x) < radius(ball)` (open ball)
or `norm(x) <= radius(ball)` (closed ball).
"""
abstract type HyperBall{T,C} <: Domain{T} end

"A ball in a fixed N-dimensional Euclidean space."
const EuclideanHyperBall{N,T,C} = HyperBall{SVector{N,T},C}

"A ball with vector elements of variable length."
const VectorHyperBall{T,C} = HyperBall{Vector{T},C}

iscompatible(x::AbstractVector, d::VectorHyperBall) = length(x) == dimension(d)

const ClosedHyperBall{T} = HyperBall{T,:closed}
const OpenHyperBall{T} = HyperBall{T,:open}

isclosedset(::ClosedHyperBall) = true
isclosedset(::OpenHyperBall) = false
isopenset(ball::HyperBall) = !isclosedset(ball)

indomain(x, ball::OpenHyperBall) = norm(x) < radius(ball)
indomain(x, ball::ClosedHyperBall) = norm(x) <= radius(ball)
approx_indomain(x, ball::OpenHyperBall, tolerance) = norm(x) < radius(ball)+tolerance
approx_indomain(x, ball::ClosedHyperBall, tolerance) = norm(x) <= radius(ball)+tolerance

# A closed ball always contains at least the origin.
isempty(ball::ClosedHyperBall) = false
isempty(ball::OpenHyperBall) = radius(ball) == 0

==(d1::HyperBall, d2::HyperBall) = isclosedset(d1)==isclosedset(d2) &&
    radius(d1)==radius(d2) && dimension(d1)==dimension(d2)

convert(::Type{SubLevelSet}, d::HyperBall{T,C}) where {T,C} =
    SubLevelSet{T,C}(norm, radius(d))
convert(::Type{SubLevelSet{T}}, d::HyperBall{S,C}) where {T,S,C} =
    SubLevelSet{T,C}(norm, radius(d))

"The unit ball."
abstract type UnitHyperBall{T,C} <: HyperBall{T,C} end

radius(::UnitHyperBall) = 1

"The unit ball with fixed dimension(s) specified by the element type."
struct StaticUnitBall{T,C} <: UnitHyperBall{T,C}
end

StaticUnitBall{T}() where {T} = StaticUnitBall{T,:closed}()

"The unit ball in a fixed N-dimensional space."
const EuclideanUnitBall{N,T,C} = StaticUnitBall{SVector{N,T},C}

EuclideanUnitBall{N}() where {N} = EuclideanUnitBall{N,Float64}()

const UnitDisk{T} = EuclideanUnitBall{2,T,:closed}
const UnitBall{T} = EuclideanUnitBall{3,T,:closed}

UnitDisk() = UnitDisk{Float64}()
UnitBall() = UnitBall{Float64}()


"The unit ball with variable dimension stored in a data field."
struct DynamicUnitBall{T,C} <: UnitHyperBall{T,C}
    dimension   ::  Int
end

DynamicUnitBall{T}(dimension::Int) where {T} =
    DynamicUnitBall{T,:closed}(dimension)

dimension(ball::DynamicUnitBall) = ball.dimension

"The unit ball with vector elements of a given dimension."
const VectorUnitBall{T,C} = DynamicUnitBall{Vector{T},C}

VectorUnitBall(dimension::Int = 3) = VectorUnitBall{Float64}(dimension)
VectorUnitDisk() = VectorUnitBall(2)

similardomain(ball::StaticUnitBall{S,C}, ::Type{T}) where {S,T,C} =
    StaticUnitBall{T,C}()
similardomain(ball::DynamicUnitBall{S,C}, ::Type{T}) where {S,T,C} =
    DynamicUnitBall{T,C}(ball.dimension)


show(io::IO, d::UnitHyperBall{T,:closed}) where {T} =
    print(io, "the $(dimension(d))-dimensional closed unit ball")
show(io::IO, d::UnitHyperBall{T,:open}) where {T} =
    print(io, "the $(dimension(d))-dimensional open unit ball")

# We choose the origin here
point_in_domain(ball::HyperBall{T}) where {T} = zero(T)
point_in_domain(ball::VectorHyperBall{T}) where {T} = zeros(T, dimension(ball))



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

isclosedset(::HyperSphere) = true
isopenset(::HyperSphere) = false

==(d1::HyperSphere, d2::HyperSphere) =
    radius(d1)==radius(d2) && dimension(d1)==dimension(d2)

convert(::Type{LevelSet}, d::HyperSphere{T}) where {T} =
    LevelSet{T}(norm, radius(d))
convert(::Type{LevelSet{T}}, d::HyperSphere) where {T} =
    LevelSet{T}(norm, radius(d))

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
VectorUnitCircle() = VectorUnitSphere(2)

similardomain(sphere::FixedUnitSphere, ::Type{T}) where {T} = FixedUnitSphere{T}()
similardomain(sphere::FlexibleUnitSphere, ::Type{T}) where {T} = FlexibleUnitSphere{T}(sphere.dimension)

show(io::IO, d::UnitHyperSphere) =
    dimension(d) == 2 ? print(io, "the unit circle") : print(io, "the $(dimension(d))-dimensional unit sphere")

point_in_domain(d::EuclideanHyperSphere{N,T}) where {N,T} =
    SVector{N,T}(ntuple( i -> i==1, N))

function point_in_domain(d::VectorHyperSphere{T}) where {T}
    p = zeros(T, dimension(d))
    p[1] = 1
    p
end

boundary(d::EuclideanUnitBall{N,T}) where {N,T} = EuclideanUnitSphere{N,T}()
boundary(d::VectorUnitBall{T}) where {T} = VectorUnitSphere{T}(dimension(d))

interior(d::EuclideanUnitBall{N,T}) where {N,T} = EuclideanUnitBall{N,T,:open}()
interior(d::VectorUnitBall{T}) where {T} = VectorUnitBall{T,:open}(dimension(d))

################
# Applications
################

"Create a cylinder with given radius and length."
cylinder(::Type{T} = Float64) where {T} = UnitDisk{T}() × UnitInterval{T}()
cylinder(radius, length) = cylinder(promote(radius, length)...)
cylinder(radius::T, length::T) where {T} = (radius .* UnitDisk{T}()) × (0..length)

"Create an ellipse curve with semi-axes lengths `a` and `b` respectively."
ellipse(a::Number, b::Number) = ellipse(promote(a,b)...)
ellipse(a::T, b::T) where {T <: Number} = scaling_map(a, b).(UnitCircle{T}())

"Create an ellipse-shaped domain with semi-axes lengths `a` and `b` respectively."
ellipse_shape(a::Number, b::Number) = ellipse_shape(promote(a,b)...)
ellipse_shape(a::T, b::T) where {T <: Number} = scaling_map(a, b).(UnitDisk{T}())


"""
The map `[cos(2πt), sin(2πt)]` from `[0,1)` to the unit circle in `ℝ^2`.
"""
struct UnitCircleMap{T} <: Map{T} end

applymap(m::UnitCircleMap{T}, t) where {T} = SVector(cos(2*T(pi)*t), sin(2*T(pi)*t))

function gradient(m::UnitCircleMap{T}, t) where {T}
    a = 2*convert(T, pi)
    SVector(-a*sin(a*t), a*cos(a*t))
end


"""
`AngleMap` is a left inverse of `UnitCircleMap`. A 2D vector `x` is projected onto
the intersection point with the unit circle of the line connecting `x` to the
origin. The angle of this point, scaled to the interval `[0,1)`, is the result.
"""
struct AngleMap{T} <: Map{SVector{2,T}}
end

function applymap(m::AngleMap{T}, x) where {T}
    twopi = 2*convert(T, pi)
    θ = atan(x[2],x[1])
    if θ < 0
        # atan2 returns an angle in (-π,π], convert to [0,2π) using periodicity.
        θ += twopi
    end
    # And divide by 2π to scale to [0,1)
    θ / twopi
end

leftinverse(m::UnitCircleMap{T}) where {T} = AngleMap{T}()
rightinverse(m::AngleMap{T}) where {T} = UnitCircleMap{T}()

canonicaldomain(d::UnitCircle{T}) where {T} = UnitInterval{T}()

fromcanonical(d::UnitCircle{T}) where {T} = UnitCircleMap{T}()
tocanonical(d::UnitCircle{T}) where {T} = AngleMap{T}()


## The complex plane

const ComplexUnitCircle{T} = FixedUnitSphere{Complex{T}}
const ComplexUnitDisk{T,C} = StaticUnitBall{Complex{T},C}

ComplexUnitCircle() = ComplexUnitCircle{Float64}()
ComplexUnitDisk() = ComplexUnitDisk{Float64}()
ComplexUnitDisk{Float64}() = ComplexUnitDisk{Float64,:closed}()

show(io::IO, d::ComplexUnitCircle{T}) where {T} =
    print(io, "the complex unit circle (T=Complex{$T})")
show(io::IO, d::ComplexUnitDisk{T,:closed}) where {T} =
    print(io, "the complex unit disk (T=Complex{$T})")
show(io::IO, d::ComplexUnitDisk{T,:open}) where {T} =
    print(io, "the complex open unit disk (T=Complex{$T})")
