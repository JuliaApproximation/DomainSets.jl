
"A `HyperRectangle` is the cartesian product of intervals."
abstract type HyperRectangle{T} <: ProductDomain{T} end

boundingbox(d::HyperRectangle) = d

"A `Cube` is a hyperrectangle with equal side lengths in each dimension."
abstract type Cube{T} <: HyperRectangle{T} end

"A cube in a fixed N-dimensional Euclidean space."
const EuclideanCube{N,T} = Cube{SVector{N,T}}

"A cube with vector elements of variable length."
const VectorCube{T} = Cube{Vector{T}}

# This will cause a warning if vectors of different length are used with cubes
iscompatiblepair(x::AbstractVector, d::VectorCube) = length(x) == dimension(d)


"The unit cube is the domain `[0,1]^d`."
abstract type UnitCube{T} <: Cube{T} end

element(d::UnitCube, i) =
    (1 <= i <= dimension(d) || throw(BoundsError); UnitInterval{numtype(d)}())

volume(d::UnitCube) = 1


"A unit cube that is specified by the element type `T`."
struct StaticUnitCube{T} <: UnitCube{T}
end

StaticUnitCube() = StaticUnitCube{Float64}()
StaticUnitCube(::Val{N}) where {N} = StaticUnitCube{SVector{N,Float64}}()

StaticUnitCube{T}(n::Int) where {T} =
    (@assert n == euclideandimension(T); StaticUnitCube{T}())
StaticUnitCube{T}(::Val{N}) where {N,T} =
    (@assert N == euclideandimension(T); StaticUnitCube{T}())

similardomain(d::StaticUnitCube, ::Type{T}) where {T<:StaticTypes} =
    StaticUnitCube{T}()
similardomain(d::StaticUnitCube, ::Type{T}) where {T} =
    DynamicUnitCube{T}(dimension(d))

elements(d::StaticUnitCube{SVector{N,T}}) where {N,T} =
    ntuple(x->UnitInterval{T}(), Val(N))

"The unit cube in a fixed N-dimensional space."
const EuclideanUnitCube{N,T} = StaticUnitCube{SVector{N,T}}

EuclideanUnitCube{N}() where {N} = EuclideanUnitCube{N,Float64}()

const UnitSquare{T} = EuclideanUnitCube{2,T}


"A unit cube whose dimension is specified by a field."
struct DynamicUnitCube{T} <: UnitCube{T}
    dimension   ::  Int

    DynamicUnitCube{T}(n::Int) where {T} = new(n)
    DynamicUnitCube{T}(n::Int) where {T<:StaticTypes} =
        (@assert n == euclideandimension(T); new(n))
end

DynamicUnitCube(n::Int) = DynamicUnitCube{Vector{Float64}}(n)

dimension(d::DynamicUnitCube) = d.dimension

elements(d::DynamicUnitCube) = map(x->UnitInterval{numtype(d)}(), 1:dimension(d))

similardomain(d::DynamicUnitCube, ::Type{T}) where {T} =
    DynamicUnitCube{T}(d.dimension)
similardomain(d::DynamicUnitCube, ::Type{T}) where {T <: StaticTypes} =
    StaticUnitCube{T}()

"The unit cube with vector elements of a given dimension."
const VectorUnitCube{T} = DynamicUnitCube{Vector{T}}

VectorUnitCube(n::Int = 3) = VectorUnitCube{Float64}(n)
VectorUnitSquare() = VectorUnitCube(2)


UnitCube(n::Int) = DynamicUnitCube(n)
UnitCube(::Val{N} = Val(3)) where {N} = EuclideanUnitCube{N}()

UnitCube{T}(n::Int) where {T <: StaticTypes} = StaticUnitCube{T}(n)
UnitCube{T}(::Val{N}) where {N,T} = StaticUnitCube{T}(Val(N))
UnitCube{T}() where {T <: StaticTypes} = StaticUnitCube{T}()
UnitCube{T}(n::Int) where {T} = DynamicUnitCube{T}(n)

UnitCube(domains::UnitInterval...) = UnitCube(domains)
UnitCube(domains::NTuple{N,UnitInterval{T}}) where {N,T} =
    UnitCube{SVector{N,T}}(domains)
UnitCube(domains::SVector{N,UnitInterval{T}}) where {N,T} =
    UnitCube{SVector{N,T}}(domains)

UnitCube{T}(domains::UnitInterval...) where {T} = UnitCube{T}(domains)
UnitCube{T}(domain::NTuple{N,<:UnitInterval}) where {N,T<:SVector{N}} =
    StaticUnitCube{T}()
UnitCube{T}(domain::SVector{N,<:UnitInterval}) where {N,T<:SVector{N}} =
    StaticUnitCube{T}()

# Constructor: careful about ambiguities with FixedInterval arguments below
ProductDomain(domains::UnitInterval...) = UnitCube(domains...)
ProductDomain(domains::NTuple{N,<:UnitInterval}) where {N} = UnitCube(domains)
ProductDomain(domains::SVector{N,<:UnitInterval}) where {N} = UnitCube(domains)
ProductDomain(domains::AbstractVector{<:UnitInterval{T}}) where {T} =
    VectorUnitCube{T}(length(domains))
ProductDomain{T}(domains::UnitInterval...) where {N,S,T<:SVector{N,S}} = UnitCube{T}(domains...)
ProductDomain{T}(domains::NTuple{N,<:UnitInterval}) where {N,S,T<:SVector{N,S}} = UnitCube{T}(domains)
ProductDomain{T}(domains::SVector{N,<:UnitInterval}) where {N,S,T<:SVector{N,S}} = UnitCube{T}(domains)
ProductDomain{T}(domains::AbstractVector{<:UnitInterval{T}}) where {T} = VectorUnitCube{T}(domains)


"An N-dimensional rectangle is the cartesian product of closed intervals."
struct Rectangle{T} <: HyperRectangle{T}
    a   ::  T
    b   ::  T

    function Rectangle{T}(a::S,b::S) where {S,T}
        @assert length(a)==length(b)
        new(a,b)
    end
end

dimension(d::Rectangle) = length(d.a)
element(d::Rectangle, i) = ClosedInterval(d.a[i],d.b[i])
elements(d::Rectangle) = map(ClosedInterval, d.a, d.b)

Rectangle(a, b) = Rectangle(promote(a,b)...)
Rectangle(a::T, b::T) where {T} = Rectangle{T}(a, b)
Rectangle(a::NTuple{N,T}, b::NTuple{N,T}) where {N,T} =
    Rectangle(SVector{N,T}(a), SVector{N,T}(b))

Rectangle(domains::Tuple) = Rectangle(domains...)
Rectangle(domains::ClosedInterval...) = Rectangle(promote_domains(domains)...)
Rectangle(domains::ClosedInterval{T}...) where {T} =
    Rectangle(map(infimum, domains), map(supremum, domains))
Rectangle(domains::AbstractVector{<:ClosedInterval}) =
    Rectangle(map(infimum, domains), map(supremum, domains))
Rectangle(domains::Domain...) =
    error("The Rectangle constructor expects two points or a list of intervals (closed).")

Rectangle{T}(domains::Tuple) where {T} = Rectangle{T}(domains...)
Rectangle{T}(domains::ClosedInterval...) where {T} =
	Rectangle{T}(map(infimum, domains), map(supremum, domains))
Rectangle{T}(domains::AbstractVector{<:ClosedInterval}) where {T} =
    Rectangle{T}(map(infimum, domains), map(supremum, domains))
Rectangle{T}(domains::Domain...) where {T} =
    error("The Rectangle constructor expects two points or a list of intervals (closed).")


ProductDomain(domains::ClosedInterval...) = Rectangle(domains...)
ProductDomain(domains::AbstractVector{<:ClosedInterval}) =
    Rectangle(map(infimum, domains), map(supremum, domains))
ProductDomain{T}(domains::ClosedInterval...) where {T} = Rectangle{T}(domains...)
ProductDomain{T}(domains::AbstractVector{<:ClosedInterval}) where {T} =
    Rectangle{T}(map(infimum, domains), map(supremum, domains))


"The N-fold cartesian product of a fixed interval."
struct FixedIntervalProduct{D,N,T} <: Cube{SVector{N,T}}
end

element(d::FixedIntervalProduct{D}, i) where {D} =
    (1 <= i <= dimension(d) || throw(BoundsError); D())
elements(d::FixedIntervalProduct{D,N}) where {D,N} =
    ntuple(x->D(), Val(N))

volume(d::FixedIntervalProduct{D,N}) where {D,N} = volume(D())^N

FixedIntervalProduct(domains::NTuple{N,D}) where {N,D <: FixedInterval} =
	FixedIntervalProduct{D,N,eltype(D)}()

const ChebyshevProductDomain{N,T} = FixedIntervalProduct{ChebyshevInterval{T},N,T}
ChebyshevProductDomain(::Val{N}) where {N} = ChebyshevProductDomain{N}()
ChebyshevProductDomain{N}() where {N} = ChebyshevProductDomain{N,Float64}()

ProductDomain(domains::D...) where {D <: FixedInterval} =
	FixedIntervalProduct(domains)
ProductDomain(domains::NTuple{N,D}) where {N,D <: FixedInterval} =
	FixedIntervalProduct(domains)
ProductDomain{T}(domains::D...) where {N,S,T<:SVector{N,S},D <: FixedInterval} =
	FixedIntervalProduct(convert.(Ref(Domain{S}), domains))
ProductDomain{T}(domains::NTuple{N,D}) where {N,S,T<:SVector{N,S},D <: FixedInterval} =
	FixedIntervalProduct(convert.(Ref(Domain{S}), domains))
