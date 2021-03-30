
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
abstract type AbstractIdentityMap{T} <: Map{T} end

applymap(map::AbstractIdentityMap, x) = x

inv(m::AbstractIdentityMap) = m

islinear(::AbstractIdentityMap) = true
isreal(::AbstractIdentityMap{T}) where {T} = eltype(T) <: Real

isidentity(::AbstractIdentityMap) = true
isidentity(m::Map{T}) where {T} = m == IdentityMap{T}()

dimension(m::AbstractIdentityMap{T}) where {T<:Number} = 1
dimension(m::AbstractIdentityMap{T}) where {N,T<:SVector{N}} = N

size(m::AbstractIdentityMap) = (dimension(m), dimension(m))

matrix(m::AbstractIdentityMap) = identitymatrix(m)
vector(m::AbstractIdentityMap) = zerovector(m)

jacobian(m::AbstractIdentityMap) = ConstantMap(matrix(m))
jacobian(m::AbstractIdentityMap, x) = matrix(m)

jacdet(m::AbstractIdentityMap, x) = 1

mapcompose(m1::AbstractIdentityMap) = m1
mapcompose(m1::AbstractIdentityMap, maps...) = mapcompose(maps...)
mapcompose2(m1, m2::AbstractIdentityMap, maps...) = mapcompose(m1, maps...)

"The identity map for variables of type `T`."
struct IdentityMap{T} <: AbstractIdentityMap{T}
end

IdentityMap() = IdentityMap{Float64}()

convert(::Type{Map{T}}, ::IdentityMap) where {T} = IdentityMap{T}()
convert(::Type{IdentityMap{T}}, ::IdentityMap) where {T} = IdentityMap{T}()

==(::IdentityMap, ::IdentityMap) = true


"Identity map with flexible size determined by a dimension field."
struct FlexibleIdentityMap{T} <: AbstractIdentityMap{T}
    dimension   ::  Int
end
const VectorIdentityMap{T} = FlexibleIdentityMap{Vector{T}}

VectorIdentityMap(dimension::Int) = VectorIdentityMap{Float64}(dimension)

dimension(m::FlexibleIdentityMap) = m.dimension

convert(::Type{Map{T}}, m::FlexibleIdentityMap) where {T} = FlexibleIdentityMap{T}(m.dimension)

==(m1::FlexibleIdentityMap, m2::FlexibleIdentityMap) = m1.dimension == m2.dimension


"The supertype of constant maps from `T` to `U`."
abstract type AbstractConstantMap{T,U} <: TypedMap{T,U} end

applymap(m::AbstractConstantMap, x) = constant(m)

isconstant(m::AbstractMap) = false
isconstant(m::AbstractConstantMap) = true

isreal(m::AbstractConstantMap) = isreal(constant(m))

dimension(m::AbstractConstantMap) = length(constant(m))
size(m::AbstractConstantMap{<:Number}) = (dimension(m), 1)
size(m::AbstractConstantMap) = (dimension(m), dimension(m))

matrix(m::AbstractConstantMap) = zeromatrix(m)
vector(m::AbstractConstantMap) = constant(m)

jacobian(m::AbstractConstantMap{T}) where {T} = ConstantMap{T}(matrix(m))
jacobian(m::AbstractConstantMap, x) = matrix(m)

jacdet(::AbstractConstantMap, x) = 0

==(m1::AbstractConstantMap, m2::AbstractConstantMap) = constant(m1)==constant(m2)

"The zero map `f(x) = 0`."
struct ZeroMap{T,U} <: AbstractConstantMap{T,U}
end

ZeroMap{T}() where {T} = ZeroMap{T,T}()

constant(m::ZeroMap{T,U}) where {T,U} = zero(T)

convert(::Type{Map{T}}, ::ZeroMap{S,U}) where {T,S,U} = ZeroMap{T,U}()

"The unity map `f(x) = 1`."
struct UnityMap{T,U} <: AbstractConstantMap{T,U}
end

UnityMap{T}() where {T} = UnityMap{T,T}()

constant(m::UnityMap{T,U}) where {T,U} = one(U)

convert(::Type{Map{T}}, ::UnityMap{S,U}) where {T,S,U} = UnityMap{T,U}()


"The constant map `f(x) = c`."
struct ConstantMap{T,U} <: AbstractConstantMap{T,U}
    c   ::  U
end

ConstantMap{T}(c::U) where {T,U} = ConstantMap{T,U}(c)
ConstantMap(c::T) where {T} = ConstantMap{T}(c)

constant(m::ConstantMap) = m.c

convert(::Type{Map{T}}, m::ConstantMap{S,U}) where {T,S,U} = ConstantMap{T,U}(m.c)


"A generic map defined by a function object."
struct GenericFunctionMap{T,F} <: Map{T}
    fun     ::  F
end

GenericFunctionMap{T}(fun) where {T} = GenericFunctionMap{T,typeof(fun)}(fun)

convert(::Type{Map{T}}, m::GenericFunctionMap{S,F}) where {S,F,T} = GenericFunctionMap{T}(m.fun)
