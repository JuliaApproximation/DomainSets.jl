
"The composition of several maps."
struct Composition{T,MAPS} <: CompositeLazyMap{T}
    maps    ::  MAPS
end

function Composition(maps...)
    T = domaintype(maps[1])
    Composition{T}(maps...)
end

Composition{T}(maps...) where {T} = Composition{T,typeof(maps)}(maps)

# TODO: make proper conversion
convert(::Type{Map{T}}, m::Composition{T}) where {T} = m
convert(::Type{Map{T}}, m::Composition) where {T} = Composition{T}(m.maps...)

applymap(m::Composition, x) = applymap_rec(x, m.maps...)
applymap_rec(x) = x
applymap_rec(x, map1, maps...) = applymap_rec(map1(x), maps...)

for op in (:inv, :leftinv, :rightinv)
    @eval $op(cmap::Composition) = Composition(reverse(map($op, elements(cmap)))...)
end

mapresulttype(m::Composition) = mapresulttype(m.maps[end])

mapcompose(m) = m
mapcompose(maps...) = Composition(maps...)

(∘)(map1::Map, map2::Map) = mapcompose(map2, map1)
(∘)(map1::Composition, map2::Map) = mapcompose(map2, elements(map1)...)
(∘)(map1::Map, map2::Composition) = mapcompose(elements(map2)..., map1)
(∘)(map1::Composition, map2::Composition) = mapcompose(elements(map2)..., elements(map1)...)


==(m1::Composition, m2::Composition) = all(map(isequal, elements(m1), elements(m2)))

## Lazy multiplication

"The lazy multiplication of one or more maps."
struct MulMap{T,MAPS} <: CompositeLazyMap{T}
    maps    ::  MAPS
end

MulMap(maps::Map{T}...) where {T} = MulMap{T}(maps...)
MulMap{T}(maps::Map{T}...) where {T} = MulMap{T,typeof(maps)}(maps)
MulMap{T}(maps...) where {T} = _mulmap(T, convert.(Map{T}, maps)...)
_mulmap(::Type{T}, maps...) where {T} = MulMap{T,typeof(maps)}(maps)

convert(::Type{Map{T}}, m::MulMap{T}) where {T} = m
convert(::Type{Map{T}}, m::MulMap) where {T} = MulMap{T}(m.maps...)

applymap(m::MulMap, x) = reduce(*, applymap.(elements(m), Ref(x)))

mapmultiply(m) = m
mapmultiply(maps...) = MulMap(maps...)


## Lazy sum

"The lazy sum of one or more maps."
struct SumMap{T,MAPS} <: CompositeLazyMap{T}
    maps    ::  MAPS
end

SumMap(maps::Map{T}...) where {T} = SumMap{T}(maps...)
SumMap{T}(maps::Map{T}...) where {T} = SumMap{T,typeof(maps)}(maps)
SumMap{T}(maps...) where {T} = _summap(T, convert.(Map{T}, maps)...)
_summap(::Type{T}, maps...) where {T} = SumMap{T,typeof(maps)}(maps)

applymap(m::SumMap, x) = reduce(+, applymap.(elements(m), Ref(x)))

mapsum(m) = m
mapsum(maps...) = SumMap(maps...)


# Define the jacobian of a composite map
jacobian(m::Composition) = composite_jacobian(elements(m)...)
composite_jacobian(map1) = jacobian(map1)
composite_jacobian(map1, map2) = mapmultiply(jacobian(map1) ∘ map2, jacobian(map2))
function composite_jacobian(map1, map2, maps...)
    rest = Composition(map2, maps...)
    f1 = jacobian(map1) ∘ rest
    f2 = composite_jacobian(map2, maps...)
    mapmultiply(f1, f2)
end

jacobian(m::MulMap) = mul_jacobian(elements(m)...)
mul_jacobian(map1) = jacobian(map1)
mul_jacobian(map1, map2) = mapsum(mapmultiply(jacobian(map1), map2), mapmultiply(map1, jacobian(map2)))
function mul_jacobian(map1, map2, maps...)
    rest = mapmultiply(map2, maps...)
    mul_jacobian(map1, rest)
end

jacobian(m::SumMap) = sum_jacobian(elements(m)...)
sum_jacobian(map1) = jacobian(map1)
sum_jacobian(maps...) = SumMap(map(jacobian, maps)...)
