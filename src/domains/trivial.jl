
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

show(io::IO, d::EmptySpace) = print(io, "{} (empty domain)")

isempty(d::EmptySpace) = true

isopenset(d::EmptySpace) = true
isclosedset(d::EmptySpace) = true

boundary(d::EmptySpace) = d
interior(d::EmptySpace) = d
closure(d::EmptySpace) = d
boundingbox(d::EmptySpace) = d

distance_to(d::EmptySpace, x) = convert(prectype(d), Inf)

# Arithmetic operations

issubset1(d1::EmptySpace, d2) = true
issubset2(d1, d2::EmptySpace) = isempty(d1)

# setdiffdomain2(d1, d2::EmptySpace) = d2

map_domain(map::Map{T}, d::EmptySpace{T}) where {T} = d
mapped_domain(map::Map, d::EmptySpace) = EmptySpace{codomaintype(map)}()

==(d1::EmptySpace, d2::EmptySpace) = true
isequal1(d1::EmptySpace, d2) = isempty(d2)
isequal2(d1, d2::EmptySpace) = isempty(d1)
hash(d::EmptySpace, h::UInt) = hash("EmptySpace", h)


"The full space of elements of type `T`."
struct FullSpace{T} <: Domain{T} end

const AnyFullSpace = FullSpace{Any}

FullSpace() = FullSpace{Float64}()
FullSpace(d) = FullSpace{eltype(d)}()

fullspace(d) = fullspace(eltype(d))
fullspace(::Type{T}) where {T} = FullSpace{T}()

isfullspace(d::FullSpace) = true
isfullspace(d::Domain) = false

similardomain(::FullSpace, ::Type{T}) where {T} = FullSpace{T}()

euclideanspace(n::Val{N}) where {N} = euclideanspace(n, Float64)
euclideanspace(::Val{N}, ::Type{T}) where {N,T} = FullSpace{SVector{N,T}}()

indomain(x::T, d::FullSpace{T}) where {T} = true
approx_indomain(x, d::FullSpace, tolerance) = in(x, d)

show(io::IO, d::FullSpace) = print(io, "{x} (full space)")

# We choose the origin as a point in the full space
point_in_domain(d::FullSpace) = zero(eltype(d))

isempty(::FullSpace) = false

isopenset(d::FullSpace) = true
isclosedset(d::FullSpace) = true

boundary(d::FullSpace{T}) where {T} = EmptySpace{T}()
interior(d::FullSpace) = d
closure(d::FullSpace) = d

distance_to(d::FullSpace, x) = zero(prectype(d))

# Arithmetic operations

issubset1(d1::FullSpace, d2) = isfullspace(d2)
issubset2(d1, d2::FullSpace) = true

map_domain(m::AbstractAffineMap{T}, d::FullSpace{T}) where {T} = d

==(d1::FullSpace, d2::FullSpace) = true
isequal1(d1::FullSpace, d2) = isfullspace(d2)
isequal2(d1, d2::FullSpace) = isfullspace(d1)
hash(d::FullSpace, h::UInt) = hash("FullSpace", h)


convert(::Type{Domain}, ::Type{T}) where T = FullSpace{T}()
convert(::Type{Domain{S}}, ::Type{T}) where {T,S} = convert(Domain{S}, convert(Domain, T))

infimum(d::FullSpace{T}) where {T} = typemin(T)
supremum(d::FullSpace{T}) where {T} = typemax(T)
