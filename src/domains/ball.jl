
# The type hierarchy is as follows:
# abstract Ball
# |-> abstract UnitBall: radius is 1
#     |-> StaticUnitBall: dimension is part of type
#     |-> DynamicUnitBall: dimension is specified by int field
# |-> GenericBall: stores center and radius. Here, the dimension can be
#     obtained from the center, so there is only one type.
# There are aliases for SVector{N,T} and Vector{T}:
#   EuclideanBall, VectorBall,
#   EuclideanUnitBall (StaticUnitBall), VectorUnitBall (DynamicUnitBall).

"""
    Ball{T,C} <: Domain{T}

Abstract supertype for volumes of elements satisfying `norm(x-center(ball)) < radius(ball)`
(open ball) or `norm(x-center(ball)) <= radius(ball)` (closed ball).
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
isopenset(d::Ball) = !isclosedset(d)

"""
    Ball(radius = 1[, center])
    Ball{T,C=:closed}(radius = 1[, center])

Return a concrete subtype of `Ball` which represents a ball with the given
radius and center, the given eltype `T`, and which is open or closed.

The default center is the origin. In case both radius and center are omitted,
a subtype of `UnitBall` is returned, whose concrete type depends on `T`.

A ball represents a volume. For the boundary of a ball, see [`Sphere()`](@ref).
"""
Ball(; radius=One(), center=Origin()) = Ball(radius, center)
Ball(radius; center=Origin()) = Ball(radius, center)
Ball(::One, ::Origin) = UnitBall()
Ball(radius, center) = GenericBall(radius, center)
Ball(::One, center) = GenericBall(1, center)
Ball(radius, center::Origin) = GenericBall(radius)

Ball{T}(args...; options...) where T = Ball{T,:closed}(args...; options...)

Ball{T,C}(; radius = One(), center = Origin()) where {T,C} = Ball{T,C}(radius, center)
Ball{T,C}(radius; center=Origin()) where {T,C} = Ball{T,C}(radius, center)
Ball{T,C}(radius, center) where {T,C} = GenericBall{T,C}(radius, center)
Ball{T,C}(::One, ::Origin) where {T,C} = UnitBall{T,C}()
Ball{T,C}(radius::One, center) where {T,C} = GenericBall{T,C}(1, center)
Ball{T,C}(radius, center::Origin) where {T,C} = GenericBall{T,C}(radius)


indomain(x, d::OpenBall) = norm(x-center(d)) < radius(d)
indomain(x, d::ClosedBall) = norm(x-center(d)) <= radius(d)
approx_indomain(x, d::OpenBall, tolerance) = norm(x-center(d)) < radius(d)+tolerance
approx_indomain(x, d::ClosedBall, tolerance) = norm(x-center(d)) <= radius(d)+tolerance

# A closed ball always contains at least its center point.
isempty(d::ClosedBall) = false
isempty(d::OpenBall) = radius(d) == 0

isequaldomain(d1::Ball, d2::Ball) = isclosedset(d1)==isclosedset(d2) &&
    radius(d1)==radius(d2) && center(d1)==center(d2)
domainhash(d::Ball, h::UInt) = hashrec("Ball", isclosedset(d), radius(d), center(d), h)

function issubset_domain(d1::Ball, d2::Ball)
    if dimension(d1) == dimension(d2)
        if center(d1) == center(d2)
            if radius(d1) < radius(d2)
                true
            elseif radius(d1) == radius(d2)
                isclosedset(d2) || isopenset(d1)
            else
                false
            end
        else
            dist = norm(center(d1)-center(d2))
            if dist + radius(d1) < radius(d2)
                true
            else
                false
            end
        end
    else
        false
    end
end

normal(d::Ball, x) = normal(boundary(d), x)

distance_to(d::Ball, x) = x ∈ d ? zero(prectype(d)) : norm(x-center(d))-radius(d)

# We choose the center of the ball here. Concrete types should implement 'center'
choice(d::Ball) = center(d)


abstract type UnitBall{T,C} <: Ball{T,C} end

radius(::UnitBall) = 1
center(d::UnitBall{T}) where {T<:StaticTypes} = zero(T)

"""
    UnitBall([dim::Int])
    UnitBall{T}([dim::Int])
    UnitBall{T,C=:closed}([dim::Int])

The open or closed volume of all points of type `T` with norm smaller than (or equal to) 1.
"""
UnitBall(::Val{N} = Val(3)) where {N} = EuclideanUnitBall{N}()
UnitBall(dim::Int) = DynamicUnitBall(dim)

UnitBall{T}(dim::Int) where {T <: StaticTypes} = StaticUnitBall{T}(dim)
UnitBall{T}(::Val{N}) where {N,T} = StaticUnitBall{T}(Val(N))
UnitBall{T}() where {T <: StaticTypes} = StaticUnitBall{T}()
UnitBall{T}(dim::Int) where {T} = DynamicUnitBall{T}(dim)

UnitBall{T,C}(dim::Int) where {T <: StaticTypes,C} = StaticUnitBall{T,C}(dim)
UnitBall{T,C}(::Val{N}) where {N,T,C} = StaticUnitBall{T,C}(Val(N))
UnitBall{T,C}() where {T <: StaticTypes,C} = StaticUnitBall{T,C}()
UnitBall{T,C}(dim::Int) where {T,C} = DynamicUnitBall{T,C}(dim)

const ClosedUnitBall{T} = UnitBall{T,:closed}
const OpenUnitBall{T} = UnitBall{T,:open}

indomain(x, d::OpenUnitBall) = norm(x) < 1
indomain(x, d::ClosedUnitBall) = norm(x) <= 1
approx_indomain(x, d::OpenUnitBall, tolerance) = norm(x) < 1+tolerance
approx_indomain(x, d::ClosedUnitBall, tolerance) = norm(x) <= 1+tolerance



isequaldomain(d1::UnitBall, d2::UnitBall) = isclosedset(d1)==isclosedset(d2) &&
    dimension(d1) == dimension(d2)

issubset_domain(d1::UnitBall, d2::UnitBall) =
    dimension(d1) == dimension(d2) && (isclosedset(d2) || isopenset(d1))

convert(::Type{SublevelSet}, d::UnitBall{T,C}) where {T,C} =
    SublevelSet{T,C}(norm, radius(d))
convert(::Type{SublevelSet{T}}, d::UnitBall{S,C}) where {T,S,C} =
    SublevelSet{T,C}(norm, radius(d))


"""
    StaticUnitBall()
    StaticUnitBall{T}()
    StaticUnitBall{T,C=:closed}()

The open or closed unit ball with static dimension determined by the element type.
"""
struct StaticUnitBall{T,C} <: UnitBall{T,C}
end

StaticUnitBall() = StaticUnitBall{Float64}()
StaticUnitBall(::Val{N}) where {N} = StaticUnitBall{SVector{N,Float64}}()

StaticUnitBall{T}() where {T} = StaticUnitBall{T,:closed}()

StaticUnitBall{T}(dim::Int) where {T} =
    (@assert dim == euclideandimension(T); StaticUnitBall{T}())
StaticUnitBall{T}(::Val{N}) where {N,T} =
    (@assert N == euclideandimension(T); StaticUnitBall{T}())

StaticUnitBall{T,C}(dim::Int) where {T,C} =
    (@assert dim == euclideandimension(T); StaticUnitBall{T,C}())
StaticUnitBall{T,C}(::Val{N}) where {N,T,C} =
    (@assert N == euclideandimension(T); StaticUnitBall{T,C}())


"The unit ball in a fixed N-dimensional space."
const EuclideanUnitBall{N,T,C} = StaticUnitBall{SVector{N,T},C}

EuclideanUnitBall{N}() where {N} = EuclideanUnitBall{N,Float64}()

similardomain(d::StaticUnitBall{S,C}, ::Type{T}) where {S,C,T<:StaticTypes} =
    StaticUnitBall{T,C}()
similardomain(d::StaticUnitBall{S,C}, ::Type{T}) where {S,C,T} =
    DynamicUnitBall{T,C}(dimension(d))

"""
    UnitDisk()
    UnitDisk{T=Float64}()

The closed unit disk in 2D, with element type `SVector{2,T}`.
"""
const UnitDisk{T} = UnitBall{SVector{2,T},:closed}
UnitDisk() = UnitDisk{Float64}()

canonicaldomain(::Parameterization, d::UnitDisk{T}) where {T} = UnitSquare{T}()
mapfrom_canonical(::Parameterization, d::UnitDisk{T}) where {T} = UnitDiskMap{T}()

const Disk{T} = Ball{SVector{2,T},:closed}
Disk() = Disk{Float64}()
Disk(radius::Number) = Disk(float(radius))
Disk(radius::AbstractFloat) = Disk{typeof(radius)}(radius)
Disk(radius,center) = Disk{float(promote_type(typeof(radius),eltype(center)))}(radius, center)

## A StaticUnitBall{<:Number} equals the interval [-1,1]  (open or closed)
convert(::Type{Interval}, d::StaticUnitBall{T,:closed}) where {T <: Number} =
    ChebyshevInterval{T}()
convert(::Type{Interval}, d::StaticUnitBall{T,:open}) where {T <: Number} =
    OpenInterval{T}(-1, 1)

canonicaldomain(::Equal, d::StaticUnitBall{T}) where {T<:Number} = convert(Interval, d)


"""
    DynamicUnitBall(dim::Int)
    DynamicUnitBall{T}(dim::Int)
    DynamicUnitBall{T,C=:closed}(dim::Int)

The open or closed unit ball with variable dimension. Typically the element
type is a `Vector{T}` and `dim` specifies the length of the vectors.
"""
struct DynamicUnitBall{T,C} <: UnitBall{T,C}
    dimension   ::  Int

    DynamicUnitBall{T,C}(dim::Int) where {T,C} = new(dim)
    DynamicUnitBall{T,C}(dim::Int) where {T<:StaticTypes,C} =
        (@assert dim == euclideandimension(T); new(dim))
end

DynamicUnitBall(dim::Int) = DynamicUnitBall{Vector{Float64}}(dim)
DynamicUnitBall{T}(dim::Int) where {T} = DynamicUnitBall{T,:closed}(dim)

dimension(d::DynamicUnitBall) = d.dimension

center(d::DynamicUnitBall{T}) where {T} = zeros(eltype(T), dimension(d))

"The unit ball with vector elements of a given dimension."
const VectorUnitBall{T,C} = DynamicUnitBall{Vector{T},C}

VectorUnitBall(dim::Int = 3) = VectorUnitBall{Float64}(dim)
VectorUnitDisk() = VectorUnitBall(2)

similardomain(d::DynamicUnitBall{S,C}, ::Type{T}) where {S,C,T} =
    DynamicUnitBall{T,C}(dimension(d))
similardomain(d::DynamicUnitBall{S,C}, ::Type{T}) where {S,C,T<:StaticTypes} =
    StaticUnitBall{T,C}()

canonicaldomain(::Parameterization, d::VectorUnitBall{T}) where {T} =
    dimension(d) == 2 ? UnitSquare{T}() : d
mapfrom_canonical(::Parameterization, d::VectorUnitBall{T}) where {T} =
    dimension(d) == 2 ? UnitDiskMap{T}() : identitymap(d)


"A `GenericBall` is a ball with a given radius and center."
struct GenericBall{T,C,S} <: Ball{T,C}
    radius  ::  S
    center  ::  T
end

GenericBall() = GenericBall(1.0)
GenericBall(radius::S) where {S<:Number} = GenericBall(radius, zero(SVector{3,float(S)}))
GenericBall(radius::S, center::T) where {S<:Number,U,T<:Vector{U}} =
    GenericBall{Vector{promote_type(S,U)}}(radius, center)
GenericBall(radius::S, center::T) where {S<:Number,N,U,T<:SVector{N,U}} =
    GenericBall{SVector{N,promote_type(S,U)}}(radius, center)
GenericBall(radius::S, center::T) where {S<:Number,T<:AbstractVector} =
    GenericBall{Vector{promote_type(S,eltype(T))}}(radius, center)
GenericBall(radius::S, center::T) where {S<:Number,T<:Number} =
    GenericBall{promote_type(S,T)}(radius, center)
GenericBall{T}(radius) where {T} = GenericBall{T,:closed}(radius)
GenericBall{T}(radius, center) where {T} = GenericBall{T,:closed}(radius, center)
GenericBall{T,C}(radius) where {T,C} = GenericBall{T,C}(radius, zero(T))
GenericBall{T,C}(radius, center) where {T,C} = GenericBall{T,C,eltype(T)}(radius, center)

GenericBall(radius::S, center::Point) where {S<:Number} =
    GenericBall(radius, pointval(center))
GenericBall{T,C}(radius, center::Point) where {T,C} =
    GenericBall{T,C}(radius, pointval(center))

radius(d::GenericBall) = d.radius
center(d::GenericBall) = d.center

similardomain(d::GenericBall{S,C}, ::Type{T}) where {S,T,C} =
    GenericBall{T,C}(d.radius, d.center)

dimension(d::GenericBall) = length(d.center)

canonicaldomain(d::GenericBall{T,C}) where {T,C} = UnitBall{T,C}(dimension(d))
mapfrom_canonical(d::GenericBall) = AffineMap(radius(d), center(d))

boundingbox(d::GenericBall) =
    map_boundingbox(boundingbox(canonicaldomain(d)), mapfrom_canonical(d))

# Preserve the `Ball` type under affine maps which preserve shape
map_domain(m::GenericAffineMap{T,S}, d::Ball{U,C}) where {T<:AbstractVector,S<:Number,U<:AbstractVector,C} =
    Ball{U,C}(m.A * radius(d), m.A * center(d) + m.b)
map_domain(m::AffineMap{T}, d::Ball{U,C}) where {T<:Number,U<:Number,C} =
    Ball{promote_type(U,T),C}(m.A * radius(d), m.A * center(d) + m.b)
map_domain(m::GenericLinearMap{T,S}, d::Ball{U,C}) where {S<:Number,T<:AbstractVector,U<:AbstractVector,C} =
    Ball{promote_type(T,U),C}(m.A * radius(d), m.A * center(d))
map_domain(m::LinearMap{T}, d::Ball{U,C}) where {T<:Number,U<:Number,C} =
    Ball{promote_type(T,U),C}(m.A * radius(d), m.A * center(d))
map_domain(m::LinearMap{T}, d::Ball{U,C}) where {T<:Number,U,C} =
    Ball{U,C}(m.A * radius(d), m.A * center(d))
map_domain(m::Translation{T}, d::Ball{S,C}) where {T<:AbstractVector,S<:AbstractVector,C} =
    Ball{S,C}(radius(d), center(d) + m.b)
map_domain(m::Translation{T}, d::Ball{S,C}) where {T<:Number,S<:Number,C} =
    Ball{S,C}(radius(d), center(d) + m.b)


show(io::IO, d::EuclideanUnitBall{3,Float64,:closed}) = print(io, "UnitBall()")
show(io::IO, d::EuclideanUnitBall{N,Float64,:closed}) where {N} = print(io, "UnitBall(Val($(N)))")
show(io::IO, d::EuclideanUnitBall{3,Float64,:open}) = print(io, "UnitBall()  (open)")
show(io::IO, d::EuclideanUnitBall{N,Float64,:open}) where {N} = print(io, "UnitBall(Val($(N)))  (open)")
Display.object_parentheses(d::EuclideanUnitBall{N,Float64,:open}) where {N} = true
show(io::IO, d::EuclideanUnitBall{2,Float64,:closed}) = print(io, "UnitDisk()")
show(io::IO, d::UnitDisk{T}) where {T} = print(io, "UnitDisk{$(T)}()")
show(io::IO, d::UnitBall{T}) where {T<:Number} = print(io, "UnitBall{$(T)}()")
show(io::IO, d::UnitBall{T,:open}) where {T<:Number} = print(io, "UnitBall{$(T)}()  (open)")
Display.object_parentheses(d::UnitBall{T,:open}) where {T<:Number} = true
show(io::IO, d::VectorUnitBall{Float64,:closed}) = print(io, "UnitBall($(dimension(d)))")
show(io::IO, d::VectorUnitBall{Float64,:open}) = print(io, "UnitBall($(dimension(d)))  (open)")
show(io::IO, d::Ball) = print(io, "Ball($(radius(d)), $(center(d)))")


# The type hierarchy of spheres parallels that of Ball above:
# abstract Sphere
# |-> abstract UnitSphere: radius is 1
#     |-> StaticUnitSphere: dimension is part of type
#     |-> DynamicUnitSphere: dimension is specified by int field
# |-> GenericSphere: stores center and radius. Here, the dimension can be
#     obtained from the center, so there is only one type.
# There are aliases for SVector{N,T} and Vector{T}.

"Supertype of spherical domains for which elements satisfy `norm(x) == radius(sphere)`."
abstract type Sphere{T} <: Domain{T} end

indomain(x, d::Sphere) = norm(x-center(d)) == radius(d)
approx_indomain(x, d::Sphere, tolerance) =
    radius(d)-tolerance <= norm(x-center(d)) <= radius(d)+tolerance

"""
    Sphere(radius = 1[, center])
    Sphere{T}(radius = 1[, center])

Return a concrete subtype of `Sphere` which represents a sphere with the given
radius and center, and the given eltype `T`.

The default center is the origin. In case both radius and center are omitted,
a subtype of `UnitSphere` is returned, whose concrete type depends on `T`.

A sphere represents the boundary of a ball. For the volume, see [`Ball()`](@ref).
"""
Sphere(; radius=One(), center=Origin()) = Sphere(radius, center)
Sphere(radius; center=Origin()) = Sphere(radius, center)
Sphere(::One, ::Origin) = UnitSphere()
Sphere(radius, center) = GenericSphere(radius, center)
Sphere(::One, center) = GenericSphere(1, center)
Sphere(radius, center::Origin) = GenericSphere(radius)

Sphere{T}(; radius = One(), center = Origin()) where T = Sphere{T}(radius, center)
Sphere{T}(radius; center=Origin()) where T = Sphere{T}(radius, center)
Sphere{T}(radius, center) where T = GenericSphere{T}(radius, center)
Sphere{T}(::One, ::Origin) where T = UnitSphere{T}()
Sphere{T}(radius::One, center) where T = GenericSphere{T}(1, center)
Sphere{T}(radius, center::Origin) where T = GenericSphere{T}(radius)

isempty(::Sphere) = false

isclosedset(::Sphere) = true
isopenset(::Sphere) = false

normal(d::Sphere, x) = (x-center(d))/norm(x-center(d))

distance_to(d::Sphere, x) = abs(norm(x-center(d))-radius(d))

choice(d::Sphere) = center(d) + unitvector(d, 1)

isequaldomain(d1::Sphere, d2::Sphere) =
    radius(d1)==radius(d2) && center(d1)==center(d2)
domainhash(d::Sphere, h::UInt) = hashrec("Sphere", radius(d), center(d), h)

"A hypersphere in a fixed N-dimensional Euclidean space."
const EuclideanSphere{N,T} = Sphere{SVector{N,T}}

"A sphere with vector elements of variable length."
const VectorSphere{T} = Sphere{Vector{T}}


"The unit sphere."
abstract type UnitSphere{T} <: Sphere{T} end

radius(::UnitSphere) = 1
center(d::UnitSphere{T}) where {T<:StaticTypes} = zero(T)

normal(d::UnitSphere, x) = x/norm(x)

indomain(x, d::UnitSphere) = norm(x) == 1
approx_indomain(x, d::UnitSphere, tolerance) =
    1-tolerance <= norm(x) <= 1+tolerance


"""
    UnitSphere([dim::Int])
    UnitSphere{T}([dim::Int])
    UnitSphere{T,C=:closed}([dim::Int])

The set of all points of type `T` with norm 1.
"""
UnitSphere(::Val{N} = Val(3)) where {N} = EuclideanUnitSphere{N}()
UnitSphere(dim::Int) = DynamicUnitSphere(dim)

UnitSphere{T}(dim::Int) where {T <: StaticTypes} = StaticUnitSphere{T}(dim)
UnitSphere{T}(::Val{N}) where {N,T} = StaticUnitSphere{T}(Val(N))
UnitSphere{T}() where {T <: StaticTypes} = StaticUnitSphere{T}()
UnitSphere{T}(dim::Int) where {T} = DynamicUnitSphere{T}(dim)

isequaldomain(d1::UnitSphere, d2::UnitSphere) = dimension(d1)==dimension(d2)

issubset1(d1::UnitSphere, d2::UnitBall) =
    dimension(d1) == dimension(d2) && isclosedset(d2)

convert(::Type{LevelSet}, d::UnitSphere{T}) where {T} = LevelSet{T}(norm, radius(d))
convert(::Type{LevelSet{T}}, d::UnitSphere) where {T} = LevelSet{T}(norm, radius(d))


"The unit sphere with fixed dimension(s) specified by the element type."
struct StaticUnitSphere{T} <: UnitSphere{T}
end

StaticUnitSphere() = StaticUnitSphere{Float64}()
StaticUnitSphere(::Val{N}) where {N} = StaticUnitSphere{SVector{N,Float64}}()

StaticUnitSphere{T}(dim::Int) where {T} =
    (@assert dim == euclideandimension(T); StaticUnitSphere{T}())
StaticUnitSphere{T}(::Val{N}) where {N,T} =
    (@assert N == euclideandimension(T); StaticUnitSphere{T}())

similardomain(d::StaticUnitSphere, ::Type{T}) where {T<:StaticTypes} =
    StaticUnitSphere{T}()
similardomain(d::StaticUnitSphere, ::Type{T}) where {T} =
    DynamicUnitSphere{T}(dimension(d))


"The unit sphere in a fixed N-dimensional Euclidean space."
const EuclideanUnitSphere{N,T} = StaticUnitSphere{SVector{N,T}}

EuclideanUnitSphere{N}() where {N} = EuclideanUnitSphere{N,Float64}()


"""
    UnitCircle()
    UnitCircle{T=Float64}()

The unit circle in 2D, with element type `SVector{2,T}`.
"""
const UnitCircle{T} = UnitSphere{SVector{2,T}}
UnitCircle() = UnitCircle{Float64}()

canonicaldomain(::Parameterization, d::UnitCircle{T}) where {T} = UnitInterval{T}()
mapfrom_canonical(::Parameterization, d::UnitCircle{T}) where {T} = UnitCircleMap{T}()


"The unit sphere with variable dimension."
struct DynamicUnitSphere{T} <: UnitSphere{T}
    dimension   ::  Int

    DynamicUnitSphere{T}(dim::Int) where {T} = new(dim)
    DynamicUnitSphere{T}(dim::Int) where {T<:StaticTypes} =
        (@assert dim == euclideandimension(T); new(dim))
end

DynamicUnitSphere(dim::Int) = DynamicUnitSphere{Vector{Float64}}(dim)

dimension(d::DynamicUnitSphere) = d.dimension

center(d::DynamicUnitSphere{T}) where {T} = zeros(eltype(T), dimension(d))

"The unit sphere with vector elements of a given dimension."
const VectorUnitSphere{T} = DynamicUnitSphere{Vector{T}}

VectorUnitSphere(dim::Int = 3) = VectorUnitSphere{Float64}(dim)
VectorUnitCircle() = VectorUnitSphere(2)

similardomain(d::DynamicUnitSphere, ::Type{T}) where {T} =
    DynamicUnitSphere{T}(d.dimension)
similardomain(d::DynamicUnitSphere, ::Type{T}) where {T <: StaticTypes} =
    StaticUnitSphere{T}()

canonicaldomain(::Parameterization, d::VectorUnitSphere{T}) where {T} =
    dimension(d) == 2 ? UnitInterval{T}() : d
mapfrom_canonical(::Parameterization, d::VectorUnitSphere{T}) where {T} =
    dimension(d) == 2 ? UnitCircleMap{T}() : identitymap(d)


"A `GenericSphere` is a sphere with a given radius and center."
struct GenericSphere{T,S} <: Sphere{T}
    radius  ::  S
    center  ::  T
end

GenericSphere() = GenericSphere(1.0)
GenericSphere(radius::S) where {S<:Number} =
    GenericSphere(radius, zero(SVector{3,float(S)}))
GenericSphere(radius::S, center::T) where {S<:Number,U,T<:Vector{U}} =
    GenericSphere{Vector{promote_type(S,U)}}(radius, center)
GenericSphere(radius::S, center::T) where {S<:Number,N,U,T<:SVector{N,U}} =
    GenericSphere{SVector{N,promote_type(S,U)}}(radius, center)
GenericSphere(radius::S, center::T) where {S<:Number,T<:AbstractVector} =
    GenericSphere{Vector{promote_type(S,eltype(T))}}(radius, center)
GenericSphere{T}(radius) where {T} = GenericSphere{T}(radius, zero(T))
GenericSphere{T}(radius, center) where {T} =
    GenericSphere{T,eltype(T)}(radius, center)

GenericSphere(radius::S, center::Point) where {S<:Number} =
    GenericSphere(radius, pointval(center))
GenericSphere{T}(radius, center::Point) where {T} =
    GenericSphere{T}(radius, pointval(center))

radius(d::GenericSphere) = d.radius
center(d::GenericSphere) = d.center

dimension(d::GenericSphere) = length(d.center)

canonicaldomain(d::GenericSphere{T}) where {T} = UnitSphere{T}(dimension(d))
mapfrom_canonical(d::GenericSphere) = AffineMap(radius(d), center(d))

boundingbox(d::GenericSphere) =
    map_boundingbox(boundingbox(canonicaldomain(d)), mapfrom_canonical(d))

# Preserve the `Sphere` type under affine maps which preserve shape
map_domain(m::GenericAffineMap{T,S}, d::Sphere{U}) where {T<:AbstractVector,S<:Number,U<:AbstractVector} =
    Sphere{U}(m.A * radius(d), m.A * center(d) + m.b)
map_domain(m::AffineMap{T}, d::Sphere{U}) where {T<:Number,U<:Number} =
    Sphere{promote_type(U,T)}(m.A * radius(d), m.A * center(d) + m.b)
map_domain(m::GenericLinearMap{T,S}, d::Sphere{U}) where {S<:Number,T<:AbstractVector,U<:AbstractVector} =
    Sphere{promote_type(T,U)}(m.A * radius(d), m.A * center(d))
map_domain(m::LinearMap{T}, d::Sphere{U}) where {T<:Number,U<:Number} =
    Sphere{promote_type(T,U)}(m.A * radius(d), m.A * center(d))
map_domain(m::LinearMap{T}, d::Sphere{U}) where {T<:Number,U} =
    Sphere{U}(m.A * radius(d), m.A * center(d))
map_domain(m::Translation{T}, d::Sphere{S}) where {T<:AbstractVector,S<:AbstractVector} =
    Sphere{S}(radius(d), center(d) + m.b)
map_domain(m::Translation{T}, d::Sphere{S}) where {T<:Number,S<:Number} =
    Sphere{S}(radius(d), center(d) + m.b)

show(io::IO, d::EuclideanUnitSphere{3,Float64}) = print(io, "UnitSphere()")
show(io::IO, d::EuclideanUnitSphere{N,Float64}) where {N} = print(io, "UnitSphere(Val($(N)))")
show(io::IO, d::EuclideanUnitSphere{2,Float64}) = print(io, "UnitCircle()")
show(io::IO, d::UnitCircle{T}) where {T} = print(io, "UnitCircle{$(T)}()")
show(io::IO, d::VectorUnitSphere{Float64}) = print(io, "UnitSphere($(dimension(d)))")
show(io::IO, d::UnitSphere{T}) where {T<:Number} = print(io, "UnitSphere{$(T)}()")
show(io::IO, d::Sphere) = print(io, "Sphere($(radius(d)), $(center(d)))")


boundary(d::UnitBall{T}) where {T} = UnitSphere{T}(dimension(d))
boundary(d::Ball{T}) where {T} = Sphere{T}(radius(d), center(d))

interior(d::StaticUnitBall{T}) where {T} = StaticUnitBall{T,:open}()
interior(d::DynamicUnitBall{T}) where {T} = DynamicUnitBall{T,:open}(dimension(d))
interior(d::Ball{T}) where {T} = Ball{T,:open}(radius(d), center(d))

closure(d::StaticUnitBall{T}) where {T} = StaticUnitBall{T,:closed}()
closure(d::DynamicUnitBall{T}) where {T} = DynamicUnitBall{T,:closed}(dimension(d))
closure(d::Ball{T}) where {T} = Ball{T,:closed}(radius(d), center(d))

boundingbox(d::UnitBall{T}) where {T<:Number} = ChebyshevInterval{T}()
boundingbox(d::UnitBall{<:StaticVector{N,T}}) where {N,T} =
    ChebyshevProductDomain{N,T}()
boundingbox(d::UnitBall{T}) where {T} =
    Rectangle{T}(-ones(eltype(T), dimension(d)), ones(eltype(T), dimension(d)))

boundingbox(d::UnitSphere{T}) where {T<:Number} = ChebyshevInterval{T}()
boundingbox(d::UnitSphere{<:StaticVector{N,T}}) where {N,T} =
    ChebyshevProductDomain{N,T}()
boundingbox(d::UnitSphere{T}) where {T} =
    Rectangle{T}(-ones(eltype(T), dimension(d)), ones(eltype(T), dimension(d)))



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
