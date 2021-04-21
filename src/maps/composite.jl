
"The composition of several maps."
struct Composition{T,MAPS} <: CompositeLazyMap{T}
    maps    ::  MAPS
end

Composition(map1, maps...) = Composition{domaintype(map1)}(map1, maps...)
Composition{T}(maps...) where {T} = Composition{T,typeof(maps)}(maps)

# TODO: make proper conversion
similarmap(m::Composition, ::Type{T}) where {T} = Composition{T}(m.maps...)

codomaintype(m::Composition) = codomaintype(m.maps[end])

# Maps are applied in the order that they appear in m.maps
applymap(m::Composition, x) = applymap_rec(x, m.maps...)
applymap_rec(x) = x
applymap_rec(x, map1, maps...) = applymap_rec(map1(x), maps...)

size(m::Composition) = (size(m.maps[end])[1], size(m.maps[1])[2])

function jacobian(m::Composition, x)
    f, fd = backpropagate(x, reverse(components(m))...)
    fd
end
backpropagate(x, m1) = (m1(x), jacobian(m1, x))
function backpropagate(x, m2, ms...)
    f, fd = backpropagate(x, ms...)
    m2(f), jacobian(m2, f) * fd
end

for op in (:inv, :leftinverse, :rightinverse)
    @eval $op(cmap::Composition) = Composition(reverse(map($op, components(cmap)))...)
end

inverse(m::Composition, x) = inverse_rec(x, reverse(components(m))...)
inverse_rec(x) = x
inverse_rec(x, map1, maps...) = inverse_rec(leftinverse(map1, x), maps...)

leftinverse(m::Composition, x) = leftinverse_rec(x, reverse(components(m))...)
leftinverse_rec(x) = x
leftinverse_rec(x, map1, maps...) = leftinverse_rec(leftinverse(map1, x), maps...)

rightinverse(m::Composition, x) = rightinverse_rec(x, reverse(components(m))...)
rightinverse_rec(x) = x
rightinverse_rec(x, map1, maps...) = rightinverse_rec(rightinverse(map1, x), maps...)

compose_map() = ()
compose_map(m) = m
compose_map(m1, m2) = compose_map1(m1, m2)
compose_map1(m1, m2) = compose_map2(m1, m2)
compose_map2(m1, m2) = Composition(m1, m2)

compose_map(m1, m2, maps...) = compose_map(compose_map(m1, m2), maps...)

compose_map(m1::Composition, m2::Composition) =
    Composition(components(m1)..., components(m2)...)
compose_map1(m1::Composition, m2) = Composition(components(m1)..., m2)
compose_map2(m1, m2::Composition) = Composition(m1, components(m2)...)

# Arguments to ∘ should be reversed before passing on to mapcompose
(∘)(map1::Map, map2::Map) = compose_map(map2, map1)


==(m1::Composition, m2::Composition) =
    ncomponents(m1) == ncomponents(m2) && all(map(isequal, components(m1), components(m2)))

Display.combinationsymbol(m::Composition) = Display.Symbol('∘')
Display.displaystencil(m::Composition) =
    composite_displaystencil(m; reversecomponents=true)
show(io::IO, mime::MIME"text/plain", m::Composition) = composite_show(io, mime, m)

## Lazy multiplication

"The lazy multiplication of one or more maps."
struct MulMap{T,MAPS} <: CompositeLazyMap{T}
    maps    ::  MAPS
end

MulMap(maps::Map{T}...) where {T} = MulMap{T,typeof(maps)}(maps)
MulMap{T}(maps::Map{T}...) where {T} = MulMap{T,typeof(maps)}(maps)
MulMap{T}(maps...) where {T} = _mulmap(T, convert.(Map{T}, maps)...)
_mulmap(::Type{T}, maps...) where {T} = MulMap{T,typeof(maps)}(maps)

similarmap(m::MulMap, ::Type{T}) where {T} = MulMap{T}(m.maps...)

applymap(m::MulMap, x) = reduce(*, applymap.(components(m), Ref(x)))

multiply_map() = ()
multiply_map(m) = m
multiply_map(m1, m2) = multiply_map1(m1, m2)
multiply_map1(m1, m2) = multiply_map2(m1, m2)
multiply_map2(m1, m2) = MulMap(m1, m2)

multiply_map(m1, m2, maps...) = multiply_map(multiply_map(m1, m2), maps...)

multiply_map(m1::MulMap, m2::MulMap) =
    MulMap(components(m1)..., components(m2)...)
multiply_map1(m1::MulMap, m2) = MulMap(components(m1)..., m2)
multiply_map2(m1, m2::MulMap) = MulMap(m1, components(m2)...)

function Display.displaystencil(m::MulMap)
    A = Any[]
    list = components(m)
    push!(A, "x -> ")
    push!(A, Display.SymbolObject(list[1]))
    push!(A, "(x)")
    for i in 2:length(list)
        push!(A, " * ")
        push!(A, Display.SymbolObject(list[i]))
        push!(A, "(x)")
    end
    A
end
show(io::IO, mime::MIME"text/plain", m::MulMap) = composite_show(io, mime, m)


## Lazy sum

"The lazy sum of one or more maps."
struct SumMap{T,MAPS} <: CompositeLazyMap{T}
    maps    ::  MAPS
end

SumMap(maps::Map{T}...) where {T} = SumMap{T}(maps...)
SumMap{T}(maps::Map{T}...) where {T} = SumMap{T,typeof(maps)}(maps)
SumMap{T}(maps...) where {T} = _summap(T, convert.(Map{T}, maps)...)
_summap(::Type{T}, maps...) where {T} = SumMap{T,typeof(maps)}(maps)

similarmap(m::SumMap, ::Type{T}) where {T} = SumMap{T}(m.maps...)

applymap(m::SumMap, x) = reduce(+, applymap.(components(m), Ref(x)))

sum_map() = ()
sum_map(m) = m
sum_map(m1, m2) = sum_map1(m1, m2)
sum_map1(m1, m2) = sum_map2(m1, m2)
sum_map2(m1, m2) = SumMap(m1, m2)

sum_map(m1, m2, maps...) = sum_map(sum_map(m1, m2), maps...)

sum_map(m1::SumMap, m2::SumMap) =
    SumMap(components(m1)..., components(m2)...)
sum_map1(m1::SumMap, m2) = SumMap(components(m1)..., m2)
sum_map2(m1, m2::SumMap) = SumMap(m1, components(m2)...)

function Display.displaystencil(m::SumMap)
    A = Any[]
    list = components(m)
    push!(A, "x -> ")
    push!(A, Display.SymbolObject(list[1]))
    push!(A, "(x)")
    for i in 2:length(list)
        push!(A, " + ")
        push!(A, Display.SymbolObject(list[i]))
        push!(A, "(x)")
    end
    A
end
show(io::IO, mime::MIME"text/plain", m::SumMap) = composite_show(io, mime, m)


# Define the jacobian of a composite map
jacobian(m::Composition) = composite_jacobian(reverse(components(m))...)
composite_jacobian(map1) = jacobian(map1)
composite_jacobian(map1, map2) = multiply_map(jacobian(map1) ∘ map2, jacobian(map2))
function composite_jacobian(map1, map2, maps...)
    rest = Composition(map2, maps...)
    f1 = jacobian(map1) ∘ rest
    f2 = composite_jacobian(map2, maps...)
    multiply_map(f1, f2)
end

jacobian(m::MulMap) = mul_jacobian(components(m)...)
mul_jacobian() = ()
mul_jacobian(map1) = jacobian(map1)
mul_jacobian(map1, map2) = sum_map(multiply_map(jacobian(map1), map2), multiply_map(map1, jacobian(map2)))
function mul_jacobian(map1, map2, maps...)
    rest = multiply_map(map2, maps...)
    mul_jacobian(map1, rest)
end

jacobian(m::SumMap) = sum_jacobian(components(m)...)
sum_jacobian() = ()
sum_jacobian(map1) = jacobian(map1)
sum_jacobian(maps...) = sum_map(map(jacobian, maps)...)
