# productmap.jl

"""
A product map is diagonal and acts on each of the components of x separately:
y = f(x) becomes y_i = f_i(x_i)
"""
struct ProductMap{N,MAPS} <: AbstractMap
    # maps has an indexable and iterable type, for example a tuple of maps
    maps    ::  MAPS
end

ProductMap(maps) = ProductMap{length(maps),typeof(maps)}(maps)

eltype(dmap::ProductMap) = promote_type(map(eltype, dmap.maps)...)

ndims{N}(::ProductMap{N}) = N

(m::ProductMap)(x) = forward_map(m, x)

elements(map::ProductMap) = map.maps
element(map::ProductMap, range::Range) = ProductMap(map.maps[range])

tensorproduct(map1::AbstractMap, map2::AbstractMap) = ProductMap((map1,map2))
tensorproduct(map1::ProductMap, map2::AbstractMap) = ProductMap((elements(map1)...,map2))
tensorproduct(map1::AbstractMap, map2::ProductMap) = ProductMap((map1,elements(map2)...))
tensorproduct(map1::ProductMap, map2::ProductMap) = ProductMap((elements(map1)...,elements(map2)...))

forward_map{N}(dmap::ProductMap{N}, x) = SVector{N}(map(forward_map, dmap.maps, x))
inverse_map{N}(dmap::ProductMap{N}, x) = SVector{N}(map(inverse_map, dmap.maps, x))

inv(dmap::ProductMap) = ProductMap(map(inv, dmap.maps))

for op in (:is_linear, :isreal)
    @eval $op(dmap::ProductMap) = reduce(&, map($op, elements(dmap)))
end


# TODO: implement jacobian
# jacobian(map::ProductMap, x) =
