
"Return a zero vector of the same size as the codomain type of the map."
zerovector(m::Map) = zerovector(m, codomaintype(m))
zerovector(m::Map, ::Type{U}) where {U} = zero(U)
zerovector(m::Map, ::Type{SVector{N,U}}) where {N,U} = zero(SVector{N,U})
# If the output type is a vector, the map itself should store the size information.
zerovector(m::Map, ::Type{<:AbstractVector{U}}) where {U} = zeros(U, size(m,1))

"Return an identity matrix with the dimensions of the map."
identitymatrix(m::Map) = identitymatrix(m, codomaintype(m))
identitymatrix(m::Map, ::Type{T}) where {T} = one(T)
identitymatrix(m::Map, ::Type{SVector{N,T}}) where {N,T} = one(SMatrix{N,N,T})
identitymatrix(m::Map, ::Type{<:AbstractVector{T}}) where {T} = Diagonal{T}(ones(size(m,1)))

"Return a zero matrix of the same size as the map."
zeromatrix(m::Map) = zeromatrix(m, matrixtype(m))
zeromatrix(m::Map, ::Type{M}) where {M} = zero(M)
zeromatrix(m::Map, ::Type{M}) where {M <: StaticArray} = zero(M)
zeromatrix(m::Map, ::Type{M}) where {M <: AbstractArray} = zeros(M, size(m))


"Supertype of identity maps."
abstract type IdentityMap{T} <: Map{T} end

IdentityMap(n::Int) = DynamicIdentityMap(n)
IdentityMap() = StaticIdentityMap()
IdentityMap(::Val{N}) where {N} = StaticIdentityMap(Val(N))

IdentityMap{T}(n::Int) where {T} = DynamicIdentityMap{T}(n)
IdentityMap{T}(n::Int) where {T<:StaticTypes} = StaticIdentityMap{T}()
IdentityMap{T}(::Val{N}) where {N,T} = StaticIdentityMap{T}(Val(N))
IdentityMap{T}() where {T} = StaticIdentityMap{T}()

applymap(map::IdentityMap, x) = x
applymap!(y, map::IdentityMap, x) = y .= x

inv(m::IdentityMap) = m

islinear(::IdentityMap) = true
isreal(::IdentityMap{T}) where {T} = isreal(T)

isidentity(::IdentityMap) = true
isidentity(m::Map{T}) where {T} = m == StaticIdentityMap{T}()

size(m::IdentityMap) = (dimension(m), dimension(m))

matrix(m::IdentityMap) = identitymatrix(m)
vector(m::IdentityMap) = zerovector(m)

jacobian(m::IdentityMap) = ConstantMap(matrix(m))
jacobian(m::IdentityMap, x) = matrix(m)

jacdet(m::IdentityMap, x) = 1

mapcompose(m1::IdentityMap) = m1
mapcompose(m1::IdentityMap, maps...) = mapcompose(maps...)
mapcompose2(m1, m2::IdentityMap, maps...) = mapcompose(m1, maps...)

show(io::IO, m::IdentityMap{T}) where {T} = print(io, "x -> x")
Display.object_parentheses(m::IdentityMap) = true

"The identity map for variables of type `T`."
struct StaticIdentityMap{T} <: IdentityMap{T}
end

StaticIdentityMap() = StaticIdentityMap{Float64}()
StaticIdentityMap(::Val{N}) where {N} = StaticIdentityMap{SVector{N,Float64}}()

StaticIdentityMap{T}(n::Int) where {T} =
    (@assert n == euclideandimension(T); StaticIdentityMap{T}())
StaticIdentityMap{T}(::Val{N}) where {N,T} =
    (@assert N == euclideandimension(T); StaticIdentityMap{T}())

similarmap(m::StaticIdentityMap, ::Type{T}) where {T<:StaticTypes} = StaticIdentityMap{T}()
similarmap(m::StaticIdentityMap, ::Type{T}) where {T} =
    DynamicIdentityMap{T}(euclideandimension(T))

convert(::Type{StaticIdentityMap{T}}, ::StaticIdentityMap) where {T} = StaticIdentityMap{T}()

==(m1::StaticIdentityMap, m2::StaticIdentityMap) = true


"Identity map with dynamic size determined by a dimension field."
struct DynamicIdentityMap{T} <: IdentityMap{T}
    dimension   ::  Int
end

const EuclideanIdentityMap{N,T} = StaticIdentityMap{SVector{N,T}}
const VectorIdentityMap{T} = DynamicIdentityMap{Vector{T}}

DynamicIdentityMap(dimension::Int) = VectorIdentityMap(dimension)
VectorIdentityMap(dimension::Int) = VectorIdentityMap{Float64}(dimension)

dimension(m::DynamicIdentityMap) = m.dimension

similarmap(m::DynamicIdentityMap, ::Type{T}) where {T} =
    DynamicIdentityMap{T}(dimension(m))
similarmap(m::DynamicIdentityMap, ::Type{T}) where {T<:StaticTypes} =
    StaticIdentityMap{T}()

==(m1::DynamicIdentityMap, m2::DynamicIdentityMap) = dimension(m1) == dimension(m2)


"The supertype of constant maps from `T` to `U`."
abstract type ConstantMap{T,U} <: TypedMap{T,U} end

applymap(m::ConstantMap, x) = constant(m)

isconstant(m::AbstractMap) = false
isconstant(m::ConstantMap) = true

isreal(m::ConstantMap{T,U}) where {T,U} =
    isreal(T) && isreal(U) && isreal(constant(m))

dimension(m::ConstantMap) = length(constant(m))
size(m::ConstantMap{<:Number}) = (dimension(m), 1)
size(m::ConstantMap) = (dimension(m), dimension(m))

matrix(m::ConstantMap) = zeromatrix(m)
vector(m::ConstantMap) = constant(m)

jacobian(m::ConstantMap{T}) where {T} = ConstantMap{T}(matrix(m))
jacobian(m::ConstantMap, x) = matrix(m)

jacdet(::ConstantMap, x) = 0

==(m1::ConstantMap, m2::ConstantMap) = constant(m1)==constant(m2)

similarmap(m::ConstantMap, ::Type{T}) where {T} = ConstantMap{T}(constant(m))
similarmap(m::ConstantMap, ::Type{T}, ::Type{U}) where {T,U} = ConstantMap{T,U}(m.c)

ConstantMap() = ConstantMap{Float64}()
ConstantMap(c) = FixedConstantMap(c)
ConstantMap{T}() where {T} = UnityMap{T}()
ConstantMap{T}(c) where {T} = FixedConstantMap{T}(c)
ConstantMap{T,U}() where {T,U} = UnityMap{T,U}()
ConstantMap{T,U}(c) where {T,U} = FixedConstantMap{T,U}(c)

show(io::IO, m::ConstantMap{T}) where {T} = print(io, "x -> $(constant(m))")
Display.object_parentheses(m::ConstantMap) = true


"The zero map `f(x) = 0`."
struct ZeroMap{T,U} <: ConstantMap{T,U}
end
ZeroMap{T}() where {T} = ZeroMap{T,T}()
constant(m::ZeroMap{T,U}) where {T,U} = zero(U)
similarmap(m::ZeroMap{S,U}, ::Type{T}) where {T,S,U} = ZeroMap{T,U}()
similarmap(m::ZeroMap, ::Type{T}, ::Type{U}) where {T,U} = ZeroMap{T,U}()


"The unity map `f(x) = 1`."
struct UnityMap{T,U} <: ConstantMap{T,U}
end
UnityMap{T}() where {T} = UnityMap{T,T}()
constant(m::UnityMap{T,U}) where {T,U} = one(U)
similarmap(m::UnityMap{S,U}, ::Type{T}) where {T,S,U} = UnityMap{T,U}()
similarmap(m::UnityMap, ::Type{T}, ::Type{U}) where {T,U} = UnityMap{T,U}()


"The constant map `f(x) = c`."
struct FixedConstantMap{T,U} <: ConstantMap{T,U}
    c   ::  U
end
FixedConstantMap{T}(c::U) where {T,U} = FixedConstantMap{T,U}(c)
FixedConstantMap(c::T) where {T} = FixedConstantMap{T}(c)
constant(m::FixedConstantMap) = m.c
