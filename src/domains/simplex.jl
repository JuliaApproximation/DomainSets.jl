
###########################
# An n-dimensional simplex
###########################

abstract type Simplex{T,C} <: Domain{T} end

isclosedset(::Simplex{T,:closed}) where {T} = true
isclosedset(::Simplex{T,:open}) where {T} = false

isopenset(d::Simplex) = !isclosedset(d)

abstract type UnitSimplex{T,C} <: Simplex{T,C} end
const ClosedUnitSimplex{T} = UnitSimplex{T,:closed}
const OpenUnitSimplex{T} = UnitSimplex{T,:open}

UnitSimplex(n::Int) = DynamicUnitSimplex(n)
UnitSimplex(::Val{N}) where {N} = EuclideanUnitSimplex{N}()

UnitSimplex{T}(n::Int) where {T <: StaticTypes} = StaticUnitSimplex{T}(n)
UnitSimplex{T}(::Val{N}) where {N,T} = StaticUnitSimplex{T}(Val(N))
UnitSimplex{T}() where {T <: StaticTypes} = StaticUnitSimplex{T}()
UnitSimplex{T}(n::Int) where {T} = DynamicUnitSimplex{T}(n)

UnitSimplex{T,C}(n::Int) where {T <: StaticTypes,C} = StaticUnitSimplex{T,C}(n)
UnitSimplex{T,C}(::Val{N}) where {N,T,C} = StaticUnitSimplex{T,C}(Val(N))
UnitSimplex{T,C}() where {T <: StaticTypes,C} = StaticUnitSimplex{T,C}()
UnitSimplex{T,C}(n::Int) where {T,C} = DynamicUnitSimplex{T,C}(n)

insimplex_closed(x) = mapreduce( t-> t >= 0, &, x) && norm(x,1) <= 1
insimplex_open(x) = mapreduce( t-> t > 0, &, x) && norm(x,1) < 1
insimplex_closed(x, tol) = mapreduce( t-> t >= -tol, &, x) && norm(x,1) <= 1+tol
insimplex_open(x, tol) = mapreduce( t-> t > -tol, &, x) && norm(x,1) < 1+tol

indomain(x, ::ClosedUnitSimplex) = insimplex_closed(x)
indomain(x, ::OpenUnitSimplex) = insimplex_open(x)
approx_indomain(x, ::ClosedUnitSimplex, tolerance) = insimplex_closed(x, tolerance)
approx_indomain(x, ::OpenUnitSimplex, tolerance) = insimplex_open(x, tolerance)

isempty(::UnitSimplex) = false

==(d1::UnitSimplex, d2::UnitSimplex) =
    isclosedset(d1)==isclosedset(d2) && dimension(d1)==dimension(d2)

boundingbox(d::UnitSimplex{T}) where {T} = UnitCube{T}(dimension(d))



struct StaticUnitSimplex{T,C} <: UnitSimplex{T,C}
end

StaticUnitSimplex() = StaticUnitSimplex{Float64}()
StaticUnitSimplex(::Val{N}) where {N} = StaticUnitSimplex{SVector{N,Float64}}()

StaticUnitSimplex{T}() where {T} = StaticUnitSimplex{T,:closed}()

StaticUnitSimplex{T}(n::Int) where {T} =
    (@assert n == euclideandimension(T); StaticUnitSimplex{T}())
StaticUnitSimplex{T}(::Val{N}) where {N,T} =
    (@assert N == euclideandimension(T); StaticUnitSimplex{T}())

StaticUnitSimplex{T,C}(n::Int) where {T,C} =
    (@assert n == euclideandimension(T); StaticUnitSimplex{T,C}())
StaticUnitSimplex{T,C}(::Val{N}) where {N,T,C} =
    (@assert N == euclideandimension(T); StaticUnitSimplex{T,C}())

const EuclideanUnitSimplex{N,T,C} = StaticUnitSimplex{SVector{N,T},C}

EuclideanUnitSimplex{N}() where {N} = EuclideanUnitSimplex{N,Float64}()


struct DynamicUnitSimplex{T,C} <: UnitSimplex{T,C}
    dimension   ::  Int

    DynamicUnitSimplex{T,C}(n::Int) where {T,C} = new(n)
    DynamicUnitSimplex{T,C}(n::Int) where {T<:StaticTypes,C} =
        (@assert n == euclideandimension(T); new(n))
end

DynamicUnitSimplex(n::Int) = DynamicUnitSimplex{Vector{Float64}}(n)
DynamicUnitSimplex{T}(n::Int) where {T} = DynamicUnitSimplex{T,:closed}(n)

dimension(d::DynamicUnitSimplex) = d.dimension

const VectorUnitSimplex{T,C} = DynamicUnitSimplex{Vector{T},C}

VectorUnitSimplex(dimension) = VectorUnitSimplex{Float64}(dimension)

show(io::IO, d::EuclideanUnitSimplex{N,Float64,:closed}) where {N} = print(io, "UnitSimplex(Val($(N)))")
show(io::IO, d::VectorUnitSimplex{Float64,:closed}) = print(io, "UnitSimplex($(dimension(d)))")


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
point_in_domain(d::UnitSimplex) = center(d)

similardomain(d::StaticUnitSimplex{S,C}, ::Type{T}) where {S,T,C} =
    StaticUnitSimplex{T,C}()
similardomain(d::DynamicUnitSimplex{S,C}, ::Type{T}) where {S,T,C} =
    DynamicUnitSimplex{T,C}(d.dimension)
