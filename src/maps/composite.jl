
"The composition of several maps."
struct CompositionMap{T,MAPS} <: CompositeLazyMap{T}
    maps    ::  MAPS
end

function CompositionMap(maps...)
    MAPS = typeof(maps)
    T = domaintype(maps[end])
    CompositionMap{T}(maps...)
end

CompositionMap{T}(maps...) where {T} = CompositionMap{T,typeof(maps)}(maps)

# TODO: make proper conversion
convert(::Type{Map{T}}, m::CompositionMap{T}) where {T} = m
convert(::Type{Map{T}}, m::CompositionMap) where {T} = CompositionMap{T}(m.maps...)

applymap(m::CompositionMap, x) = applymap_rec(x, m.maps...)
applymap_rec(x) = x
applymap_rec(x, map1, maps...) = applymap_rec(map1(x), maps...)

for op in (:inv, :leftinv, :rightinv)
    @eval $op(cmap::CompositionMap) = CompositionMap(reverse(map($op, elements(cmap)))...)
end

(∘)(map1::Map, map2::Map) = CompositionMap(map2, map1)
(∘)(map1::CompositionMap, map2::Map) = CompositionMap(map2, elements(map1)...)
(∘)(map1::Map, map2::CompositionMap) = CompositionMap(elements(map2)..., map1)
(∘)(map1::CompositionMap, map2::CompositionMap) = CompositionMap(elements(map2)..., elements(map1)...)


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


"The lazy sum of one or more maps."
struct SumMap{T,MAPS} <: CompositeLazyMap{T}
    maps    ::  MAPS
end

SumMap(maps::Map{T}...) where {T} = SumMap{T}(maps...)
SumMap{T}(maps::Map{T}...) where {T} = SumMap{T,typeof(maps)}(maps)
SumMap{T}(maps...) where {T} = _summap(T, convert.(Map{T}, maps)...)
_summap(::Type{T}, maps...) where {T} = SumMap{T,typeof(maps)}(maps)

applymap(m::SumMap, x) = reduce(+, applymap.(elements(m), Ref(x)))

# Define the jacobian of a composite map
jacobian(m::CompositionMap) = composite_jacobian(elements(m)...)
composite_jacobian(map1) = jacobian(map1)
composite_jacobian(map1, map2) = MulMap(jacobian(map1) ∘ map2, jacobian(map2))
function composite_jacobian(map1, map2, maps...)
    rest = CompositionMap(map2, maps...)
    f1 = jacobian(map1) ∘ rest
    f2 = composite_jacobian(map2, maps...)
    MulMap(f1, f2)
end

jacobian(m::MulMap) = mul_jacobian(elements(m)...)
mul_jacobian(map1) = jacobian(map1)
mul_jacobian(map1, map2) = SumMap(MulMap(jacobian(map1), map2), MulMap(map1, jacobian(map2)))
function mul_jacobian(map1, map2, maps...)
    rest = MulMap(map2, maps...)
    mul_jacobian(map1, rest)
end

jacobian(m::SumMap) = sum_jacobian(elements(m)...)
sum_jacobian(map1) = jacobian(map1)
sum_jacobian(maps...) = SumMap(map(jacobian, maps)...)
