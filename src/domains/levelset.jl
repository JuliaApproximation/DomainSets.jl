
"A `FunctionLevelSet` is a set that derives from the levels of a function."
abstract type FunctionLevelSet{T} <: Domain{T} end

# Convenience: assume the function and the level are stored fields
levelfun(d::FunctionLevelSet) = d.f
levelfun(d::FunctionLevelSet, x) = d.f(x)
level(d::FunctionLevelSet) = d.level


"Supertype of level set domains of the form `f(x)=C`."
abstract type AbstractLevelSet{T} <: FunctionLevelSet{T} end

indomain(x, d::AbstractLevelSet) = levelfun(d, x) == level(d)

show(io::IO, d::AbstractLevelSet) =
    print(io, "level set f(x) = $(level(d)) with f = $(levelfun(d))")


"The domain defined by `f(x)=0` for a given function `f`."
struct ZeroSet{T,F} <: AbstractLevelSet{T}
    f   ::  F
end

ZeroSet(f) = ZeroSet{Float64}(f)
ZeroSet{T}(f::F) where {T,F} = ZeroSet{T,F}(f)

level(d::ZeroSet) = 0

similardomain(d::ZeroSet, ::Type{T}) where {T} = ZeroSet{T}(levelfun(d))


"The domain defined by `f(x)=C` for a given function `f` and constant `C`."
struct LevelSet{T,F,S} <: AbstractLevelSet{T}
    f       ::  F
    level   ::  S
end

LevelSet(f, level) = LevelSet{typeof(level)}(f, level)
LevelSet{T}(f::F, level::S) where {T,F,S} = LevelSet{T,F,S}(f, level)

similardomain(d::LevelSet, ::Type{T}) where {T} = LevelSet{T}(levelfun(d), level(d))

convert(::Type{LevelSet}, d::ZeroSet{T}) where {T} = LevelSet{T}(levelfun(d), level(d))


"Supertype of sublevel set domains."
abstract type AbstractSubLevelSet{T,C} <: FunctionLevelSet{T} end

indomain(x, d::AbstractSubLevelSet{T,:closed}) where {T} = levelfun(d, x) <= level(d)
indomain(x, d::AbstractSubLevelSet{T,:open}) where {T} = levelfun(d, x) < level(d)

show(io::IO, d::AbstractSubLevelSet{T,:closed}) where {T} =
    print(io, "sublevel set f(x) <= $(level(d)) with f = $(levelfun(d))")
show(io::IO, d::AbstractSubLevelSet{T,:open}) where {T} =
    print(io, "sublevel set f(x) < $(level(d)) with f = $(levelfun(d))")

"The domain where `f(x) <= 0` (or `f(x) < 0`)."
struct SubZeroSet{T,C,F} <: AbstractSubLevelSet{T,C}
    f   ::  F
end

SubZeroSet(f) = SubZeroSet{Float64}(f)
SubZeroSet{T}(f) where {T} = SubZeroSet{T,:closed}(f)
SubZeroSet{T,C}(f::F) where {T,C,F} = SubZeroSet{T,C,F}(f)

level(d::SubZeroSet) = 0

"The domain defined by `f(x) <= C` (or `f(x) < C`) for a given function `f` and constant `C`."
struct SubLevelSet{T,C,F,S} <: AbstractSubLevelSet{T,C}
    f       ::  F
    level   ::  S
end

SubLevelSet(f, level) = SubLevelSet{typeof(level)}(f, level)
SubLevelSet{T}(f, level) where {T} = SubLevelSet{T,:closed}(f, level)
SubLevelSet{T,C}(f::F, level::S) where {T,C,F,S} = SubLevelSet{T,C,F,S}(f, level)


"Supertype of superlevel set domains."
abstract type AbstractSuperLevelSet{T,C} <: FunctionLevelSet{T} end

indomain(x, d::AbstractSuperLevelSet{T,:closed}) where {T} = levelfun(d, x) >= level(d)
indomain(x, d::AbstractSuperLevelSet{T,:open}) where {T} = levelfun(d, x) > level(d)

show(io::IO, d::AbstractSuperLevelSet{T,:closed}) where {T} =
    print(io, "sublevel set f(x) >= $(level(d)) with f = $(levelfun(d))")
show(io::IO, d::AbstractSuperLevelSet{T,:open}) where {T} =
    print(io, "sublevel set f(x) > $(level(d)) with f = $(levelfun(d))")

"The domain where `f(x) >= 0` (or `f(x) > 0`)."
struct SuperZeroSet{T,C,F} <: AbstractSuperLevelSet{T,C}
    f   ::  F
end

SuperZeroSet(f) = SuperZeroSet{Float64}(f)
SuperZeroSet{T}(f) where {T} = SuperZeroSet{T,:closed}(f)
SuperZeroSet{T,C}(f::F) where {T,C,F} = SuperZeroSet{T,C,F}(f)

level(d::SuperZeroSet) = 0

"The domain defined by `f(x) >= C` (or `f(x) > C`) for a given function `f` and constant `C`."
struct SuperLevelSet{T,C,F,S} <: AbstractSuperLevelSet{T,C}
    f       ::  F
    level   ::  S
end

SuperLevelSet(f, level) = SuperLevelSet{typeof(level)}(f, level)
SuperLevelSet{T}(f, level) where {T} = SuperLevelSet{T,:closed}(f, level)
SuperLevelSet{T,C}(f::F, level::S) where {T,C,F,S} = SuperLevelSet{T,C,F,S}(f, level)


## Additional functionality

pseudolevel(d::AbstractLevelSet, epsilon) =
    _pseudolevel(d, epsilon, levelfun(d), level(d))
_pseudolevel(d::AbstractLevelSet{T}, epsilon, fun, C) where {T} =
    SubLevelSet{T,:open}(x -> norm(fun(x)-C), epsilon)

pseudolevel(d::Domain, epsilon) = pseudolevel(convert(LevelSet, d), epsilon)
