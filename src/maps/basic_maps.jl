# basic_maps.jl

"""
The identity map.
"""
struct IdentityMap{T} <: AbstractMap{T,T}
end

(m::IdentityMap{T})(x) where {T} = applymap(m, x)

applymap(map::IdentityMap{T}, x::T) where {T} = x

apply_inverse(map::IdentityMap{T}, y::T) where T = y

inv(m::IdentityMap) = m

islinear(::IdentityMap) = true



"""
The constant map `f(x) = c`.
"""
struct ConstantMap{T,S} <: AbstractMap{T,S}
    c   ::  T
end

constant(m::ConstantMap) = m.c

applymap(m::ConstantMap, x) = constant(m)
