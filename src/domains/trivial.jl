
"The empty domain with elements of type `T`."
struct EmptySpace{T} <: Domain{T}
end

const AnyEmptySpace = EmptySpace{Any}

EmptySpace() = EmptySpace{Float64}()
EmptySpace(::Type{T}) where {T} = EmptySpace{T}()

similardomain(::EmptySpace, ::Type{T}) where {T} = EmptySpace{T}()

emptyspace(d) = emptyspace(eltype(d))
emptyspace(::Type{T}) where {T} = EmptySpace{T}()

indomain(x::T, d::EmptySpace{T}) where {T} = false

approx_indomain(x, d::EmptySpace, tolerance) = in(x, d)

isempty(d::EmptySpace) = true

isopenset(d::EmptySpace) = true
isclosedset(d::EmptySpace) = true

boundary(d::EmptySpace) = d
interior(d::EmptySpace) = d
closure(d::EmptySpace) = d

# Arithmetic operations

union(d1::EmptySpace, d2::EmptySpace) = d1
union(d1::Domain, d2::EmptySpace) = d1
union(d1::EmptySpace, d2::Domain) = d2

intersect(d1::EmptySpace, d2::EmptySpace) = d1
intersect(d1::Domain, d2::EmptySpace) = d2
intersect(d1::EmptySpace, d2::Domain) = d1

setdiff(d1::EmptySpace, d2::EmptySpace) = d1
setdiff(d1::EmptySpace, d2::Domain) = d1
setdiff(d1::Domain, d2::EmptySpace) = d1

map_domain(map::Map{T}, d::EmptySpace{T}) where {T} = d
mapped_domain(map::Map, d::EmptySpace) = EmptySpace{codomaintype(map)}()

==(::EmptySpace, ::EmptySpace) = true

show(io::IO, d::EmptySpace) = print(io, "{} (empty domain)")


"The full space of elements of type `T`."
struct FullSpace{T} <: Domain{T} end

const AnyFullSpace = FullSpace{Any}

FullSpace() = FullSpace{Float64}()
FullSpace(d) = FullSpace{eltype(d)}()

fullspace(d) = fullspace(eltype(d))
fullspace(::Type{T}) where {T} = FullSpace{T}()


similardomain(::FullSpace, ::Type{T}) where {T} = FullSpace{T}()

euclideanspace(n::Val{N}) where {N} = euclideanspace(n, Float64)
euclideanspace(::Val{N}, ::Type{T}) where {N,T} = FullSpace{SVector{N,T}}()

indomain(x::T, d::FullSpace{T}) where {T} = true

approx_indomain(x, d::FullSpace, tolerance) = in(x, d)

# We choose the origin as a point in the full space
point_in_domain(d::FullSpace) = zero(eltype(d))

isempty(::FullSpace) = false

isopenset(d::FullSpace) = true
isclosedset(d::FullSpace) = true

boundary(d::FullSpace{T}) where {T} = EmptySpace{T}()
interior(d::FullSpace) = d
closure(d::FullSpace) = d

# Arithmetic operations

union(d1::FullSpace, d2::FullSpace) = d1
union(d1::Domain, d2::FullSpace) = d2
union(d1::FullSpace, d2::Domain) = d1

intersect(d1::FullSpace, d2::FullSpace) = d1
intersect(d1::Domain, d2::FullSpace) = d1
intersect(d1::FullSpace, d2::Domain) = d2


map_domain(m::AbstractAffineMap{T}, d::FullSpace{T}) where {T} = d

show(io::IO, d::FullSpace) = print(io, "{x} (full space)")


convert(::Type{Domain}, ::Type{T}) where T = FullSpace{T}()
convert(::Type{Domain{S}}, ::Type{T}) where {T,S} = convert(Domain{S}, convert(Domain, T))

infimum(d::FullSpace{T}) where {T} = typemin(T)
supremum(d::FullSpace{T}) where {T} = typemax(T)


# Some convenient complete domains

"The set of integers of type Int (ℤ = \\BbbZ)."
const ℤ = FullSpace{Int}()
"The set of rational numbers of type Rational{Int} (ℚ = \\BbbQ)."
const ℚ = FullSpace{Rational{Int}}()
"The set of reals of type Float64 (ℝ = \\BbbR)."
const ℝ = FullSpace{Float64}()
"The complex plane with Float64 real and imaginar parts (ℂ = \\BbbC)."
const ℂ = FullSpace{Complex{Float64}}()

"The space ℝ^1"
const ℝ1 = FullSpace{SVector{1,Float64}}()
"The space ℝ^2"
const ℝ2 = FullSpace{SVector{2,Float64}}()
"The space ℝ^3"
const ℝ3 = FullSpace{SVector{3,Float64}}()
"The space ℝ^4"
const ℝ4 = FullSpace{SVector{4,Float64}}()
