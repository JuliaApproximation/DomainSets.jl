# basic_maps.jl

"""
The identity map.
"""
struct IdentityMap{T} <: AbstractMap{T,T}
end

(m::IdentityMap)(x) = applymap(m, x)

applymap(map::IdentityMap, x) = x

apply_inverse(map::IdentityMap, y) = y

inv(m::IdentityMap) = m
