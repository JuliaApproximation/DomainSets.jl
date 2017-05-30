# composite_map.jl

"""
The composition of several maps.
"""
struct CompositeMap{MAPS} <: AbstractMap
    maps    ::  MAPS
end

(m::CompositeMap)(x) = forward_map(m, x)

elements(map::CompositeMap) = map.maps
element(map::CompositeMap, i::Int) = map.maps[i]
element(map::CompositeMap, range::Range) = CompositeMap(map.maps[range])

# dest_type{T}(map::CompositeMap, ::Type{T}) = promote_type(map(m->dest_type(m,T), elements(map))...)

forward_map(map::CompositeMap, x) = forward_map_rec(x, map.maps...)

forward_map_rec(x) = x
forward_map_rec(x, map1::AbstractMap, maps::AbstractMap...) = forward_map_rec(map1*x, maps...)

inv(cmap::CompositeMap) = CompositeMap(reverse(map(inv, cmap.maps)))

for op in (:is_linear, :isreal)
    @eval $op(cmap::CompositeMap) = reduce(&, map($op, elements(cmap)))
end

# TODO: implement jacobian
# jacobian(map::CompositeMap, x) =

(*)(map1::AbstractMap, map2::AbstractMap) = CompositeMap((map2,map1))
(*)(map1::CompositeMap, map2::AbstractMap) = CompositeMap((map2, elements(map1)...))
(*)(map1::AbstractMap, map2::CompositeMap) = CompositeMap((elements(map2)..., map1))
(*)(map1::CompositeMap, map2::CompositeMap) = CompositeMap((elements(map2)...,elements(map1)...))
