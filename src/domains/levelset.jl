
"Supertype of level set domains."
abstract type AbstractLevelSet{T} <: Domain{T} end

indomain(x, d::AbstractLevelSet) = d.f(x) == level(d)

pseudolevel(d::AbstractLevelSet, epsilon) =
    _pseudolevel(d, epsilon, levelfun(d), level(d))
_pseudolevel(d::AbstractLevelSet{T}, epsilon, fun, C) where {T} =
    IndicatorFunction{T}(x -> norm(fun(x)-C) < epsilon)

show(io::IO, d::AbstractLevelSet) =
    print(io, "level set f(x) = $(level(d)) with f = $(levelfun(d))")


"The domain defined by `f(x)=0` for a given function `f`."
struct ZeroSet{T,F} <: AbstractLevelSet{T}
    f   ::  F
end

ZeroSet(f) = ZeroSet{Float64}(f)
ZeroSet{T}(f::F) where {T,F} = ZeroSet{T,F}(f)

levelfun(d::ZeroSet) = d.f
level(d::ZeroSet) = 0

similardomain(d::ZeroSet, ::Type{T}) where {T} = ZeroSet{T}(d.f)


"The domain defined by `f(x)=C` for a given function `f` and constant `C`."
struct LevelSet{T,F,S} <: AbstractLevelSet{T}
    f   ::  F
    C   ::  S
end

LevelSet(f, C) = LevelSet{typeof(C)}(f, C)
LevelSet{T}(f::F, C::S) where {T,F,S} = LevelSet{T,F,S}(f, C)

levelfun(d::LevelSet) = d.f
level(d::LevelSet) = d.C

similardomain(d::LevelSet, ::Type{T}) where {T} = LevelSet{T}(d.f, d.C)
