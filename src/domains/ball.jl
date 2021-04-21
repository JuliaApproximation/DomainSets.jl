
# The type hierarchy is as follows:
# abstract Ball
# |-> abstract UnitBall: radius is 1
#     |-> StaticUnitBall: dimension is part of type
#     |-> DynamicUnitBall: dimension is specified by int field
# There are aliases for SVector{N,T} and Vector{T}:
#   EuclideanBall, VectorBall,
#   EuclideanUnitBall (StaticUnitBall), VectorUnitBall (DynamicUnitBall).

"""
Supertype of balls for which elements satisfy `norm(x) < radius(ball)` (open ball)
or `norm(x) <= radius(ball)` (closed ball).
"""
abstract type Ball{T,C} <: Domain{T} end

"A ball in a fixed N-dimensional Euclidean space."
const EuclideanBall{N,T,C} = Ball{SVector{N,T},C}

"A ball with vector elements of variable length."
const VectorBall{T,C} = Ball{Vector{T},C}

iscompatiblepair(x::AbstractVector, d::VectorBall) = length(x) == dimension(d)

const ClosedBall{T} = Ball{T,:closed}
const OpenBall{T} = Ball{T,:open}

isclosedset(::ClosedBall) = true
isclosedset(::OpenBall) = false
isopenset(ball::Ball) = !isclosedset(ball)

indomain(x, ball::OpenBall) = norm(x) < radius(ball)
indomain(x, ball::ClosedBall) = norm(x) <= radius(ball)
approx_indomain(x, ball::OpenBall, tolerance) = norm(x) < radius(ball)+tolerance
approx_indomain(x, ball::ClosedBall, tolerance) = norm(x) <= radius(ball)+tolerance

# A closed ball always contains at least the origin.
isempty(ball::ClosedBall) = false
isempty(ball::OpenBall) = radius(ball) == 0

==(d1::Ball, d2::Ball) = isclosedset(d1)==isclosedset(d2) &&
    radius(d1)==radius(d2) && dimension(d1)==dimension(d2)

convert(::Type{SublevelSet}, d::Ball{T,C}) where {T,C} =
    SublevelSet{T,C}(norm, radius(d))
convert(::Type{SublevelSet{T}}, d::Ball{S,C}) where {T,S,C} =
    SublevelSet{T,C}(norm, radius(d))

"The unit ball."
abstract type UnitBall{T,C} <: Ball{T,C} end

radius(::UnitBall) = 1

UnitBall(n::Int) = DynamicUnitBall(n)
UnitBall(::Val{N} = Val(3)) where {N} = EuclideanUnitBall{N}()

UnitBall{T}(n::Int) where {T <: StaticTypes} = StaticUnitBall{T}(n)
UnitBall{T}(::Val{N}) where {N,T} = StaticUnitBall{T}(Val(N))
UnitBall{T}() where {T <: StaticTypes} = StaticUnitBall{T}()
UnitBall{T}(n::Int) where {T} = DynamicUnitBall{T}(n)

UnitBall{T,C}(n::Int) where {T <: StaticTypes,C} = StaticUnitBall{T,C}(n)
UnitBall{T,C}(::Val{N}) where {N,T,C} = StaticUnitBall{T,C}(Val(N))
UnitBall{T,C}() where {T <: StaticTypes,C} = StaticUnitBall{T,C}()
UnitBall{T,C}(n::Int) where {T,C} = DynamicUnitBall{T,C}(n)


"The unit ball with fixed dimension(s) specified by the element type."
struct StaticUnitBall{T,C} <: UnitBall{T,C}
end

StaticUnitBall() = StaticUnitBall{Float64}()
StaticUnitBall(::Val{N}) where {N} = StaticUnitBall{SVector{N,Float64}}()

StaticUnitBall{T}() where {T} = StaticUnitBall{T,:closed}()

StaticUnitBall{T}(n::Int) where {T} =
    (@assert n == euclideandimension(T); StaticUnitBall{T}())
StaticUnitBall{T}(::Val{N}) where {N,T} =
    (@assert N == euclideandimension(T); StaticUnitBall{T}())

StaticUnitBall{T,C}(n::Int) where {T,C} =
    (@assert n == euclideandimension(T); StaticUnitBall{T,C}())
StaticUnitBall{T,C}(::Val{N}) where {N,T,C} =
    (@assert N == euclideandimension(T); StaticUnitBall{T,C}())

"The unit ball in a fixed N-dimensional space."
const EuclideanUnitBall{N,T,C} = StaticUnitBall{SVector{N,T},C}

EuclideanUnitBall{N}() where {N} = EuclideanUnitBall{N,Float64}()

similardomain(d::StaticUnitBall{S,C}, ::Type{T}) where {S,C,T<:StaticTypes} =
    StaticUnitBall{T,C}()
similardomain(d::StaticUnitBall{S,C}, ::Type{T}) where {S,C,T} =
    DynamicUnitBall{T,C}(dimension(d))

const UnitDisk{T} = EuclideanUnitBall{2,T,:closed}
UnitDisk() = UnitDisk{Float64}()

## A StaticUnitBall{<:Number} equals the interval [-1,1]  (open or closed)
convert(::Type{Interval}, d::StaticUnitBall{T,:closed}) where {T <: Number} =
    ChebyshevInterval{T}()
convert(::Type{Interval}, d::StaticUnitBall{T,:open}) where {T <: Number} =
    OpenInterval{T}(-1, 1)

canonicaldomain(d::StaticUnitBall{T}, ::Equal) where {T<:Number} = convert(Interval, d)

# canonicaldomain(d::StaticUnitBall{SVector{1,T}}, ::Isomorphic) where {T} =
    # convert(Interval, convert(Domain{T}, d))

"The unit ball with variable dimension stored in a data field."
struct DynamicUnitBall{T,C} <: UnitBall{T,C}
    dimension   ::  Int

    DynamicUnitBall{T,C}(n::Int) where {T,C} = new(n)
    DynamicUnitBall{T,C}(n::Int) where {T<:StaticTypes,C} =
        (@assert n == euclideandimension(T); new(n))
end

DynamicUnitBall(n::Int) = DynamicUnitBall{Vector{Float64}}(n)
DynamicUnitBall{T}(n::Int) where {T} = DynamicUnitBall{T,:closed}(n)

dimension(ball::DynamicUnitBall) = ball.dimension

"The unit ball with vector elements of a given dimension."
const VectorUnitBall{T,C} = DynamicUnitBall{Vector{T},C}

VectorUnitBall(n::Int = 3) = VectorUnitBall{Float64}(n)
VectorUnitDisk() = VectorUnitBall(2)

similardomain(d::DynamicUnitBall{S,C}, ::Type{T}) where {S,C,T} =
    DynamicUnitBall{T,C}(dimension(d))
similardomain(d::DynamicUnitBall{S,C}, ::Type{T}) where {S,C,T<:StaticTypes} =
    StaticUnitBall{T,C}()


show(io::IO, d::EuclideanUnitBall{3,Float64,:closed}) = print(io, "UnitBall()")
show(io::IO, d::EuclideanUnitBall{N,Float64,:closed}) where {N} = print(io, "UnitBall(Val($(N)))")
show(io::IO, d::EuclideanUnitBall{3,Float64,:open}) = print(io, "UnitBall()  (open)")
show(io::IO, d::EuclideanUnitBall{N,Float64,:open}) where {N} = print(io, "UnitBall(Val($(N)))  (open)")
show(io::IO, d::UnitDisk{Float64}) = print(io, "UnitDisk()")
show(io::IO, d::UnitDisk{T}) where {T} = print(io, "UnitDisk{$(T)}()")
show(io::IO, d::VectorUnitBall{Float64,:closed}) = print(io, "UnitBall($(dimension(d)))")
show(io::IO, d::VectorUnitBall{Float64,:open}) = print(io, "UnitBall($(dimension(d)))  (open)")

# We choose the origin here
point_in_domain(ball::Ball{T}) where {T} = zero(T)
point_in_domain(ball::VectorBall{T}) where {T} = zeros(T, dimension(ball))


# The type hierarchy of spheres parallels that of Ball above:
# abstract Sphere
# |-> abstract UnitSphere: radius is 1
#     |-> StaticUnitSphere: dimension is part of type
#     |-> DynamicUnitSphere: dimension is specified by int field
# There are aliases for SVector{N,T} and Vector{T}.

"Supertype of spherical domains for which elements satisfy `norm(x) == radius(sphere)`."
abstract type Sphere{T} <: Domain{T} end

indomain(x, d::Sphere) = norm(x) == radius(d)
approx_indomain(x, d::Sphere, tolerance) =
    radius(d)-tolerance <= norm(x) <= radius(d)+tolerance

isempty(::Sphere) = false

isclosedset(::Sphere) = true
isopenset(::Sphere) = false

==(d1::Sphere, d2::Sphere) =
    radius(d1)==radius(d2) && dimension(d1)==dimension(d2)

convert(::Type{LevelSet}, d::Sphere{T}) where {T} =
    LevelSet{T}(norm, radius(d))
convert(::Type{LevelSet{T}}, d::Sphere) where {T} =
    LevelSet{T}(norm, radius(d))

"A hypersphere in a fixed N-dimensional Euclidean space."
const EuclideanSphere{N,T} = Sphere{SVector{N,T}}

"A sphere with vector elements of variable length."
const VectorSphere{T} = Sphere{Vector{T}}


"The unit sphere."
abstract type UnitSphere{T} <: Sphere{T} end

radius(::UnitSphere) = 1

UnitSphere(n::Int) = DynamicUnitSphere(n)
UnitSphere(::Val{N} = Val(3)) where {N} = EuclideanUnitSphere{N}()

UnitSphere{T}(n::Int) where {T <: StaticTypes} = StaticUnitSphere{T}(n)
UnitSphere{T}(::Val{N}) where {N,T} = StaticUnitSphere{T}(Val(N))
UnitSphere{T}() where {T <: StaticTypes} = StaticUnitSphere{T}()
UnitSphere{T}(n::Int) where {T} = DynamicUnitSphere{T}(n)

issubset1(d1::UnitSphere, d2::UnitBall) = dimension(d1) == dimension(d2)

"The unit sphere with fixed dimension(s) specified by the element type."
struct StaticUnitSphere{T} <: UnitSphere{T}
end

StaticUnitSphere() = StaticUnitSphere{Float64}()
StaticUnitSphere(::Val{N}) where {N} = StaticUnitSphere{SVector{N,Float64}}()

StaticUnitSphere{T}(n::Int) where {T} =
    (@assert n == euclideandimension(T); StaticUnitSphere{T}())
StaticUnitSphere{T}(::Val{N}) where {N,T} =
    (@assert N == euclideandimension(T); StaticUnitSphere{T}())

similardomain(d::StaticUnitSphere, ::Type{T}) where {T<:StaticTypes} =
    StaticUnitSphere{T}()
similardomain(d::StaticUnitSphere, ::Type{T}) where {T} =
    DynamicUnitSphere{T}(dimension(d))


"The unit sphere in a fixed N-dimensional Euclidean space."
const EuclideanUnitSphere{N,T} = StaticUnitSphere{SVector{N,T}}

EuclideanUnitSphere{N}() where {N} = EuclideanUnitSphere{N,Float64}()


"The unit circle in 2D."
const UnitCircle{T} = EuclideanUnitSphere{2,T}

"The unit sphere with variable dimension."
struct DynamicUnitSphere{T} <: UnitSphere{T}
    dimension   ::  Int

    DynamicUnitSphere{T}(n::Int) where {T} = new(n)
    DynamicUnitSphere{T}(n::Int) where {T<:StaticTypes} =
        (@assert n == euclideandimension(T); new(n))
end

DynamicUnitSphere(n::Int) = DynamicUnitSphere{Vector{Float64}}(n)

dimension(d::DynamicUnitSphere) = d.dimension

"The unit sphere with vector elements of a given dimension."
const VectorUnitSphere{T} = DynamicUnitSphere{Vector{T}}

VectorUnitSphere(n::Int = 3) = VectorUnitSphere{Float64}(n)
VectorUnitCircle() = VectorUnitSphere(2)

similardomain(d::DynamicUnitSphere, ::Type{T}) where {T} =
    DynamicUnitSphere{T}(d.dimension)
similardomain(d::DynamicUnitSphere, ::Type{T}) where {T <: StaticTypes} =
    StaticUnitSphere{T}()

show(io::IO, d::EuclideanUnitSphere{3,Float64}) = print(io, "UnitSphere()")
show(io::IO, d::EuclideanUnitSphere{N,Float64}) where {N} = print(io, "UnitSphere(Val($(N)))")
show(io::IO, d::UnitCircle{Float64}) = print(io, "UnitCircle()")
show(io::IO, d::UnitCircle{T}) where {T} = print(io, "UnitCircle{$(T)}()")
show(io::IO, d::VectorUnitSphere{Float64}) = print(io, "UnitSphere($(dimension(d)))")

point_in_domain(d::EuclideanSphere{N,T}) where {N,T} =
    SVector{N,T}(ntuple( i -> i==1, N))

function point_in_domain(d::VectorSphere{T}) where {T}
    p = zeros(T, dimension(d))
    p[1] = 1
    p
end

boundary(d::StaticUnitBall{T}) where {T} = StaticUnitSphere{T}()
boundary(d::DynamicUnitBall{T}) where {T} = DynamicUnitSphere{T}(dimension(d))

interior(d::StaticUnitBall{T}) where {T} = StaticUnitBall{T,:open}()
interior(d::DynamicUnitBall{T}) where {T} = DynamicUnitBall{T,:open}(dimension(d))
closure(d::StaticUnitBall{T}) where {T} = StaticUnitBall{T,:closed}()
closure(d::DynamicUnitBall{T}) where {T} = DynamicUnitBall{T,:closed}(dimension(d))

boundingbox(d::UnitBall{T}) where {T<:Number} = ChebyshevInterval{T}()
boundingbox(d::UnitBall{SVector{N,T}}) where {N,T} =
    ChebyshevProductDomain{N,T}()
boundingbox(d::UnitBall{T}) where {T} =
    ProductDomain{T}((ChebyshevInterval{eltype(T)}() for i in 1:dimension(d))...)

boundingbox(d::UnitSphere{T}) where {T<:Number} = ChebyshevInterval{T}()
boundingbox(d::UnitSphere{SVector{N,T}}) where {N,T} =
    ChebyshevProductDomain{N,T}()
boundingbox(d::UnitSphere{T}) where {T} =
    ProductDomain{T}((ChebyshevInterval{eltype(T)}() for i in 1:dimension(d))...)


################
# Applications
################

"Create a cylinder with given radius and length."
cylinder(::Type{T} = Float64) where {T} = UnitDisk{T}() × UnitInterval{T}()
cylinder(radius, length) = cylinder(promote(radius, length)...)
cylinder(radius::T, length::T) where {T} = (radius .* UnitDisk{T}()) × (0..length)

"Create an ellipse curve with semi-axes lengths `a` and `b` respectively."
ellipse(a::Number, b::Number) = ellipse(promote(a,b)...)
ellipse(a::T, b::T) where {T <: Number} = LinearMap(a, b).(UnitCircle{T}())

"Create an ellipse-shaped domain with semi-axes lengths `a` and `b` respectively."
ellipse_shape(a::Number, b::Number) = ellipse_shape(promote(a,b)...)
ellipse_shape(a::T, b::T) where {T <: Number} = LinearMap(a, b).(UnitDisk{T}())


"""
The map `[cos(2πt), sin(2πt)]` from `[0,1]` to the unit circle in `ℝ^2`.
"""
struct UnitCircleMap{T} <: Map{T} end

applymap(m::UnitCircleMap{T}, t) where {T} = SVector(cos(2*T(pi)*t), sin(2*T(pi)*t))
function applymap!(y, m::UnitCircleMap{T}, t) where {T}
    y[1] = cos(2*T(pi)*t)
    y[2] = sin(2*T(pi)*t)
    y
end

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

canonicaldomain(d::UnitCircle{T}, ::Parameterization) where {T} = UnitInterval{T}()
fromcanonical(d::UnitCircle{T}, ::Parameterization) where {T} = UnitCircleMap{T}()


## The complex plane

const ComplexUnitCircle{T} = StaticUnitSphere{Complex{T}}
const ComplexUnitDisk{T,C} = StaticUnitBall{Complex{T},C}

ComplexUnitCircle() = ComplexUnitCircle{Float64}()
ComplexUnitDisk() = ComplexUnitDisk{Float64}()
ComplexUnitDisk{T}() where {T} = ComplexUnitDisk{T,:closed}()

show(io::IO, d::ComplexUnitCircle{Float64}) = print(io, "ComplexUnitCircle()")
show(io::IO, d::ComplexUnitDisk{Float64,:closed}) = print(io, "ComplexUnitDisk()")
show(io::IO, d::ComplexUnitDisk{Float64,:open}) = print(io, "ComplexUnitDisk()  (open)")
show(io::IO, d::ComplexUnitDisk{T,:closed}) where {T} = print(io, "ComplexUnitDisk{$(T)}()")
show(io::IO, d::ComplexUnitDisk{T,:open}) where {T} = print(io, "ComplexUnitDisk{$(T)}()  (open)")
