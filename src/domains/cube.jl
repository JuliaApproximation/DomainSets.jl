
"A `HyperRectangle` is the cartesian product of intervals."
abstract type AbstractHyperRectangle{T} <: ProductDomain{T} end

"A `HyperCube` is a hyperrectangle with equal side lengths in each dimension."
abstract type HyperCube{T} <: AbstractHyperRectangle{T} end

"The unit cube is the domain `[0,1]^d`."
abstract type UnitHyperCube{T} <: HyperCube{T} end

element(d::UnitHyperCube, i) =
    (1 <= i <= dimension(d) || throw(BoundsError); UnitInterval{numtype(d)}())

volume(d::UnitHyperCube) = 1

"A unit cube that is specified by the element type `T`."
struct StaticUnitCube{T} <: UnitHyperCube{T}
end

elements(d::StaticUnitCube{SVector{N,T}}) where {N,T} =
    ntuple(x->UnitInterval{T}(), Val(N))

"The unit cube in a fixed N-dimensional space."
const EuclideanUnitCube{N,T} = StaticUnitCube{SVector{N,T}}
const UnitSquare{T} = EuclideanUnitCube{2,T}
const UnitCube{T} = EuclideanUnitCube{3,T}

EuclideanUnitCube{N}() where {N} = EuclideanUnitCube{N,Float64}()
UnitSquare() = UnitSquare{Float64}()
UnitCube() = UnitCube{Float64}()

"A unit cube whose dimension is specified by a field."
struct DynamicUnitCube{T} <: UnitHyperCube{T}
    dimension   ::  Int
end

dimension(d::DynamicUnitCube) = d.dimension

"The unit cube with vector elements of a given dimension."
const VectorUnitCube{T} = DynamicUnitCube{Vector{T}}

elements(d::DynamicUnitCube) = map(x->UnitInterval{numtype(d)}(), 1:dimension(d))

UnitHyperCube(n::Int) = UnitHyperCube{Vector{T}}(n)
UnitHyperCube(domains::UnitInterval...) = UnitHyperCube(domains)
UnitHyperCube(domains::NTuple{N,UnitInterval{T}}) where {N,T} =
    UnitHyperCube{SVector{N,T}}(domains)
UnitHyperCube(domains::SVector{N,UnitInterval{T}}) where {N,T} =
    UnitHyperCube{SVector{N,T}}(domains)
UnitHyperCube{T}(n::Int) where {T} = DynamicUnitCube{T}(n)
UnitHyperCube{T}(domain::NTuple{N,<:UnitInterval}) where {N,T<:SVector{N}} =
    StaticUnitCube{T}()
UnitHyperCube{T}(domain::SVector{N,<:UnitInterval}) where {N,T<:SVector{N}} =
    StaticUnitCube{T}()

ProductDomain(domains::UnitInterval...) = UnitHyperCube(domains...)
ProductDomain(domains::NTuple{N,<:UnitInterval}) where {N} = UnitHyperCube(domains)
ProductDomain(domains::SVector{N,<:UnitInterval}) where {N} = UnitHyperCube(domains)
ProductDomain(domains::AbstractVector{<:UnitInterval{T}}) where {T} =
    UnitHyperCube{T}(length(domains))
ProductDomain{T}(domains::UnitInterval...) where {T} = UnitHyperCube{T}(domains...)
ProductDomain{T}(domains::NTuple{N,<:UnitInterval}) where {N,T} = UnitHyperCube{T}(domains)
ProductDomain{T}(domains::SVector{N,<:UnitInterval}) where {N,T} = UnitHyperCube{T}(domains)
ProductDomain{T}(domains::AbstractVector{<:UnitInterval{T}}) where {T} = UnitHyperCube{T}(domains)

"An N-dimensional hyperrectangle is the cartesian product of closed intervals."
struct HyperRectangle{T} <: AbstractHyperRectangle{T}
    a   ::  T
    b   ::  T

    function HyperRectangle{T}(a::T,b::T) where {T}
        @assert length(a)==length(b)
        new(a,b)
    end
end

dimension(d::HyperRectangle) = length(d.a)
element(d::HyperRectangle, i) = ClosedInterval(d.a[i],d.b[i])
elements(d::HyperRectangle) = map(ClosedInterval, d.a, d.b)

HyperRectangle(a, b) = HyperRectangle(promote(a,b)...)
HyperRectangle(a::T, b::T) where {T} = HyperRectangle{T}(a, b)
HyperRectangle(a::NTuple{N,T}, b::NTuple{N,T}) where {N,T} =
    HyperRectangle(SVector{N,T}(a), SVector{N,T}(b))

HyperRectangle(domains::ClosedInterval...) = HyperRectangle(promote_domains(domains)...)
HyperRectangle(domains::ClosedInterval{T}...) where {T} =
    HyperRectangle(map(infimum, domains), map(supremum, domains))
HyperRectangle(domains::AbstractVector{<:ClosedInterval}) =
    HyperRectangle(map(infimum, domains), map(supremum, domains))
HyperRectangle(domains::Domain...) =
    error("The HyperRectangle constructor expects two points or a list of intervals (closed).")

ProductDomain(domains::ClosedInterval...) = HyperRectangle(domains...)
ProductDomain(domains::AbstractVector{<:ClosedInterval}) =
    HyperRectangle(map(infimum, domains), map(supremum, domains))

_boundary(d::HyperRectangle, domains) =
	UnionDomain(ProductDomain(domains[1:i-1]..., boundary(domains[i]), domains[i+1:end]...) for i in 1:length(domains))
