# productmap.jl

"""
A product map is diagonal and acts on each of the components of x separately:
y = f(x) becomes y_i = f_i(x_i)
"""
struct ProductMap{MAPS,S <: Tuple,T <: Tuple} <: AbstractMap{S,T}
    # maps has an indexable and iterable type, for example a tuple of maps
    maps    ::  MAPS
end

function ProductMap(maps...)
    MAPS = typeof(maps)
    S = typeof(map(x->zero(domaintype(x)), maps))
    T = typeof(map(x->zero(codomaintype(x)), maps))

    ProductMap{MAPS,S,T}(maps)
end

elements(pm::ProductMap) = pm.maps

(m::ProductMap)(x) = applymap(m, x)

applymap(pm::ProductMap, x) = map(applymap, elements(pm), x)

tensorproduct(map1::AbstractMap, map2::AbstractMap) = ProductMap(map1, map2)
tensorproduct(map1::ProductMap, map2::AbstractMap) = ProductMap(elements(map1)..., map2)
tensorproduct(map1::AbstractMap, map2::ProductMap) = ProductMap(map1, elements(map2)...)
tensorproduct(map1::ProductMap, map2::ProductMap) = ProductMap(elements(map1)..., elements(map2)...)

inv(pm::ProductMap) = ProductMap(map(inv, elements(pm))...)

left_inverse(pm::ProductMap) = ProductMap(map(left_inverse, elements(pm))...)
right_inverse(pm::ProductMap) = ProductMap(map(right_inverse, elements(pm))...)
