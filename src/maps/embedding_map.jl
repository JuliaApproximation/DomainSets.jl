# embedding_map.jl

"A map between embeded spaces."
struct EmbeddingMap{T,S} <: AbstractMap{T,S}
    function EmbeddingMap{T,S}() where {T,S}
        @assert embedded(spacetype(T), spacetype(S))
        new{T,S}()
    end
end

applymap(map::EmbeddingMap{T,S}, x::T) where {T,S} = convert_space(spacetype(S), x)


"A map between isomorphic spaces."
struct IsomorphismMap{T,S} <: AbstractMap
    function IsomorphismMap{T,S}() where {T,S}
        @assert embedded(spacetype(T), spacetype(S))
        new{T,S}()
    end
end

applymap(map::IsomorphismMap{T,S}, x::T) where {T,S} = convert_space(spacetype(S), x)

inverse_map(map::IsomorphismMap{T,S}, y::S) where {T,S} = convert_space(spacetype(T), y)
