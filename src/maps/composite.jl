
"The composition of several maps."
struct CompositeMap{T,MAPS} <: Map{T}
    maps    ::  MAPS
end

function CompositeMap(maps...)
    MAPS = typeof(maps)
    T = domaintype(maps[end])
    CompositeMap{T,MAPS}(maps)
end

compose(map::Map) = map
compose(maps...) = CompositeMap(maps...)

elements(map::CompositeMap) = map.maps

applymap(map::CompositeMap, x) = applymap_rec(x, map.maps...)

applymap_rec(x) = x
applymap_rec(x, map1, maps...) = applymap_rec(map1(x), maps...)

isreal(m::CompositeMap) = all(map(isreal, elements(m)))

for op in (:inv, :leftinv, :rightinv)
    @eval $op(cmap::CompositeMap) = CompositeMap(reverse(map($op, elements(cmap)))...)
end

(∘)(map1::Map, map2::Map) = CompositeMap(map2, map1)
(∘)(map1::CompositeMap, map2::Map) = CompositeMap(map2, elements(map1)...)
(∘)(map1::Map, map2::CompositeMap) = CompositeMap(elements(map2)..., map1)
(∘)(map1::CompositeMap, map2::CompositeMap) = CompositeMap(elements(map2)..., elements(map1)...)


# struct CompositeJacobian{S,T} <: Map{S,T}
#     cmap        ::  CompositeMap
#     jacobians   ::  Vector{Map}
# end
#
# function CompositeJacobian(cmap::CompositeMap{M,S,T}) where {M,S,T}
#     jacobians = [jacobian(m) for m in elements(cmap)]
#     U = jac_type(S,T)
#     CompositeJacobian{S,U}(cmap, jacobians)
# end
#
# applymap(map::CompositeJacobian, x) = _applymap(map, x, map.jacobians, elements(map.cmap))
#
# # This function is very type-unstable
# function _applymap(map::CompositeJacobian, x, jacobians, elements)
#     r = applymap(jacobians[1], x)
#     y = applymap(elements[1], x)
#     for i in 2:length(elements)
#         r = r * applymap(jacobians[i], y)
#         y = applymap(elements[i], y)
#     end
#     r
# end
#
# jacobian(map::CompositeMap) = CompositeJacobian(map)
