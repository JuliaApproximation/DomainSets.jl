# basic_maps.jl


"""
The identity map.
"""
struct IdentityMap{T} <: AbstractMap{T,T}
end

(m::IdentityMap)(x) = applymap(m, x)

applymap(map::IdentityMap{T}, x::T) where T = x

apply_inverse(map::IdentityMap{T}, y::T) where T = y

inv(m::IdentityMap) = m

islinear(::IdentityMap) = true

jacobian(m::IdentityMap{T}) where T = UnityMap{T,jac_type(T,T)}()


"The supertype of several maps that map to a constant value."
abstract type AbstractConstantMap{S,T} <: AbstractMap{S,T}
end

applymap(m::AbstractConstantMap, x) = constant(m)

islinear(::AbstractConstantMap) = true

jacobian(m::AbstractConstantMap{S,T}) where {S,T} = ZeroMap{S,jac_type(S,T)}()


"""
The zero map `f(x) = 0`.
"""
struct ZeroMap{S,T} <: AbstractConstantMap{S,T}
end

constant(m::ZeroMap{S,T}) where {S,T} = zero(T)


"""
The unity map `f(x) = 1`.
"""
struct UnityMap{S,T} <: AbstractConstantMap{S,T}
end

constant(m::UnityMap{S,T}) where {S,T} = one(T)


"""
The constant map `f(x) = c`.
"""
struct ConstantMap{S,T} <: AbstractConstantMap{S,T}
    c   ::  T
end

ConstantMap(c::T) where T = ConstantMap{T,T}(c)

constant(m::ConstantMap) = m.c
