
###########################
# An n-dimensional simplex
###########################

abstract type Simplex{T,C} <: Domain{T} end

isclosedset(::Simplex{T,:closed}) where {T} = true
isclosedset(::Simplex{T,:open}) where {T} = false

isopenset(d::Simplex) = !isclosedset(d)

abstract type AbstractUnitSimplex{T,C} <: Simplex{T,C} end

const ClosedUnitSimplex{T} = AbstractUnitSimplex{T,:closed}
const OpenUnitSimplex{T} = AbstractUnitSimplex{T,:open}

indomain(x, ::ClosedUnitSimplex) = mapreduce( t-> t >= 0, &, x) && norm(x,1) <= 1
indomain(x, ::OpenUnitSimplex) = mapreduce( t-> t >= 0, &, x) && norm(x,1) < 1

approx_indomain(x, ::AbstractUnitSimplex, tolerance) = mapreduce( t-> t >= -tolerance, &, x) && norm(x,1) <= 1+tolerance

isempty(::AbstractUnitSimplex) = false

==(d1::AbstractUnitSimplex, d2::AbstractUnitSimplex) =
    isclosedset(d1)==isclosedset(d2) && dimension(d1)==dimension(d2)

struct FixedUnitSimplex{T,C} <: AbstractUnitSimplex{T,C}
end

FixedUnitSimplex{T}() where {T} = FixedUnitSimplex{T,:closed}()

struct FlexibleUnitSimplex{T,C} <: AbstractUnitSimplex{T,C}
    dimension   ::  Int
end

FlexibleUnitSimplex{T}(dimension) where {T} = FlexibleUnitSimplex{T,:closed}(dimension)

dimension(d::FlexibleUnitSimplex) = d.dimension

const EuclideanUnitSimplex{N,T,C} = FixedUnitSimplex{SVector{N,T},C}
const VectorUnitSimplex{T,C} = FlexibleUnitSimplex{Vector{T},C}
const UnitSimplex{N,T,C} = EuclideanUnitSimplex{N,T,C}

EuclideanUnitSimplex{N}() where {N} = EuclideanUnitSimplex{N,Float64}()
VectorUnitSimplex(dimension) = VectorUnitSimplex{Float64}(dimension)

center(d::EuclideanUnitSimplex{N,T}) where {N,T} = ones(SVector{N,T})/N
center(d::VectorUnitSimplex{T}) where {T} = ones(T, dimension(d))/dimension(d)

corners(::EuclideanUnitSimplex{N,T}) where {N,T} =
    [ SVector{N,T}(ntuple( i -> i==j, N)) for j in 1:N]
corners(d::VectorUnitSimplex{T}) where {T} =
    [ (z = zeros(T, dimension(d)); z[j]=1; z) for j in 1:dimension(d)]


# We pick the center point, because it belongs to the domain regardless of
# whether it is open or closed.
point_in_domain(d::AbstractUnitSimplex) = center(d)

convert(::Type{Domain{T}}, d::FixedUnitSimplex{S,C}) where {S,T,C} =
    FixedUnitSimplex{T,C}()
convert(::Type{Domain{T}}, d::FlexibleUnitSimplex{S,C}) where {S,T,C} =
    FlexibleUnitSimplex{T,C}(d.dimension)
