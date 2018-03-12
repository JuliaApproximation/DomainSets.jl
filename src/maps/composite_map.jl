# composite_map.jl

"""
The composition of several maps.
"""
struct CompositeMap{MAPS,T,S} <: AbstractMap{T,S}
    maps    ::  MAPS
end

function CompositeMap(maps...)
    if length(maps) == 1
        return maps[1]
    end
    MAPS = typeof(maps)
    # TODO: check all intermediate types
    S = domaintype(maps[1])
    T = rangetype(maps[end])
    CompositeMap{MAPS,T,S}(maps)
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
