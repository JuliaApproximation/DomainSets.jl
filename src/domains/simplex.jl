
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

insimplex_closed(x) = mapreduce( t-> t >= 0, &, x) && norm(x,1) <= 1
insimplex_open(x) = mapreduce( t-> t > 0, &, x) && norm(x,1) < 1
insimplex_closed(x, tol) = mapreduce( t-> t >= -tol, &, x) && norm(x,1) <= 1+tol
insimplex_open(x, tol) = mapreduce( t-> t > -tol, &, x) && norm(x,1) < 1+tol

indomain(x, ::ClosedUnitSimplex) = insimplex_closed(x)
indomain(x, ::OpenUnitSimplex) = insimplex_open(x)
approx_indomain(x, ::ClosedUnitSimplex, tolerance) = insimplex_closed(x, tolerance)
approx_indomain(x, ::OpenUnitSimplex, tolerance) = insimplex_open(x, tolerance)

isempty(::AbstractUnitSimplex) = false

==(d1::AbstractUnitSimplex, d2::AbstractUnitSimplex) =
    isclosedset(d1)==isclosedset(d2) && dimension(d1)==dimension(d2)

struct StaticUnitSimplex{T,C} <: AbstractUnitSimplex{T,C}
end

StaticUnitSimplex{T}() where {T} = StaticUnitSimplex{T,:closed}()

struct DynamicUnitSimplex{T,C} <: AbstractUnitSimplex{T,C}
    dimension   ::  Int
end

DynamicUnitSimplex{T}(dimension) where {T} = DynamicUnitSimplex{T,:closed}(dimension)

dimension(d::DynamicUnitSimplex) = d.dimension


const EuclideanUnitSimplex{N,T,C} = StaticUnitSimplex{SVector{N,T},C}
const VectorUnitSimplex{T,C} = DynamicUnitSimplex{Vector{T},C}
const UnitSimplex{N,T,C} = EuclideanUnitSimplex{N,T,C}

EuclideanUnitSimplex{N}() where {N} = EuclideanUnitSimplex{N,Float64}()
VectorUnitSimplex(dimension) = VectorUnitSimplex{Float64}(dimension)

center(d::EuclideanUnitSimplex{N,T}) where {N,T} = ones(SVector{N,T})/N
center(d::VectorUnitSimplex{T}) where {T} = ones(T, dimension(d))/dimension(d)

corners(::EuclideanUnitSimplex{N,T}) where {N,T} =
    vcat([zero(SVector{N,T})], [ SVector{N,T}(ntuple( i -> i==j, N)) for j in 1:N])
corners(d::VectorUnitSimplex{T}) where {T} =
    vcat([zeros(T,dimension(d))], [ (z = zeros(T, dimension(d)); z[j]=1; z) for j in 1:dimension(d)])

interior(d::EuclideanUnitSimplex{N,T}) where {N,T} = EuclideanUnitSimplex{N,T,:open}()
closure(d::EuclideanUnitSimplex{N,T}) where {N,T} = EuclideanUnitSimplex{N,T,:closed}()
interior(d::VectorUnitSimplex{T}) where {T} = VectorUnitSimplex{T,:open}(dimension(d))
closure(d::VectorUnitSimplex{T}) where {T} = VectorUnitSimplex{T,:closed}(dimension(d))


# We pick the center point, because it belongs to the domain regardless of
# whether it is open or closed.
point_in_domain(d::AbstractUnitSimplex) = center(d)

similardomain(d::StaticUnitSimplex{S,C}, ::Type{T}) where {S,T,C} =
    StaticUnitSimplex{T,C}()
similardomain(d::DynamicUnitSimplex{S,C}, ::Type{T}) where {S,T,C} =
    DynamicUnitSimplex{T,C}(d.dimension)
