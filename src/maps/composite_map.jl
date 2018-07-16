# composite_map.jl

"""
The composition of several maps.
"""
struct CompositeMap{MAPS,S,T} <: AbstractMap{S,T}
    maps    ::  MAPS
end

function CompositeMap(maps...)
    if length(maps) == 1
        return maps[1]
    end
    MAPS = typeof(maps)
    # TODO: check all intermediate types
    S = domaintype(maps[1])
    T = codomaintype(maps[end])
    CompositeMap{MAPS,S,T}(maps)
end


(m::CompositeMap)(x) = applymap(m, x)

elements(map::CompositeMap) = map.maps

applymap(map::CompositeMap, x) = applymap_rec(x, map.maps...)

applymap_rec(x) = x
applymap_rec(x, map1::AbstractMap, maps::AbstractMap...) = applymap_rec(map1*x, maps...)

for op in (:inv, :left_inverse, :right_inverse)
    @eval $op(cmap::CompositeMap) = CompositeMap(reverse(map($op, elements(cmap)))...)
end

(∘)(map1::AbstractMap, map2::AbstractMap) = CompositeMap(map2, map1)
(∘)(map1::CompositeMap, map2::AbstractMap) = CompositeMap(map2, elements(map1)...)
(∘)(map1::AbstractMap, map2::CompositeMap) = CompositeMap(elements(map2)..., map1)
(∘)(map1::CompositeMap, map2::CompositeMap) = CompositeMap(elements(map2)..., elements(map1)...)


struct CompositeJacobian{S,T} <: AbstractMap{S,T}
    cmap        ::  CompositeMap
    jacobians   ::  Vector{AbstractMap}
end

function CompositeJacobian(cmap::CompositeMap{M,S,T}) where {M,S,T}
    jacobians = [jacobian(m) for m in elements(cmap)]
    U = jac_type(S,T)
    CompositeJacobian{S,U}(cmap, jacobians)
end

applymap(map::CompositeJacobian, x) = _applymap(map, x, map.jacobians, elements(map.cmap))

# This function is very type-unstable
function _applymap(map::CompositeJacobian, x, jacobians, elements)
    r = applymap(jacobians[1], x)
    y = applymap(elements[1], x)
    for i in 2:length(elements)
        r = r * applymap(jacobians[i], y)
        y = applymap(elements[i], y)
    end
    r
end

jacobian(map::CompositeMap) = CompositeJacobian(map)
