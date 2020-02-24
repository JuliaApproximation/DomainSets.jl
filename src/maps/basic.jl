
"Supertype of identity maps."
abstract type AbstractIdentityMap{T} <: Map{T} end

applymap(map::AbstractIdentityMap, x) = x

inv(m::AbstractIdentityMap) = m

islinear(::AbstractIdentityMap) = true
isreal(::AbstractIdentityMap{T}) where {T} = eltype(T) <: Real

jacdet(m::AbstractIdentityMap, x) = 1

"The identity map for variables of type `T`."
struct IdentityMap{T} <: AbstractIdentityMap{T}
end

convert(::Type{Map{T}}, ::IdentityMap) where {T} = IdentityMap{T}()

"Identity map with flexible size determined by a dimension field."
struct FlexibleIdentityMap{T} <: AbstractIdentityMap{T}
    dimension   ::  Int
end
const VectorIdentityMap{T} = FlexibleIdentityMap{Vector{T}}

VectorIdentityMap(dimension::Int) = VectorIdentityMap{Float64}(dimension)

dimension(m::FlexibleIdentityMap) = m.dimension

convert(::Type{Map{T}}, m::FlexibleIdentityMap) where {T} = FlexibleIdentityMap{T}(m.dimension)

identitymatrix(::Type{SVector{N,T}}) where {N,T} = one(SMatrix{N,N,T})
identitymatrix(::Type{T}) where {T} = one(T)
identitymatrix(::Type{Vector{T}}, dimension::Int) where {T} = Matrix{T}(I, dimension, dimension)

jacobian(m::IdentityMap{T}) where {T} = ConstantMap(identitymatrix(T))
jacobian(m::IdentityMap{T}, x) where {T} = identitymatrix(T)

jacobian(m::FlexibleIdentityMap{T}) where {T} = ConstantMap(identitymatrix(T, dimension(m)))
jacobian(m::FlexibleIdentityMap{T}, x) where {T} = identitymatrix(T, dimension(m))


"The supertype of constant maps from `T` to `U`."
abstract type AbstractConstantMap{T,U} <: TypedMap{T,U} end

applymap(m::AbstractConstantMap, x) = constant(m)

islinear(::AbstractConstantMap) = true
isreal(m::AbstractConstantMap) = isreal(constant(m))

jacobian(m::AbstractConstantMap{T}) where {T} = ConstantMap{T}(zero(constant(m)))
jacobian(m::AbstractConstantMap, x) = zero(constant(m))

jacdet(::AbstractConstantMap, x) = 0


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
