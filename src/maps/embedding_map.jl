# embedding_map.jl

"A map between embedded spaces."
struct EmbeddingMap{S,T} <: AbstractMap{S,T}
    function EmbeddingMap{S,T}() where {S,T}
        @assert embedded(spacetype(S), spacetype(T))
        new{S,T}()
    end
end

embedding_map(::Type{S}, ::Type{T}) where {S,T} = EmbeddingMap{S,T}()
embedding_map(::Type{T}, ::Type{T}) where {T} = IdentityMap{T}()

(m::EmbeddingMap)(x) = applymap(m, x)

applymap(map::EmbeddingMap{S,T}, x::S) where {S,T} = convert_space(spacetype(T), x)


"A restriction map from a space to an embedded space."
struct RestrictionMap{S,T} <: AbstractMap{S,T}
    function RestrictionMap{S,T}() where {S,T}
        @assert embedded(spacetype(T), spacetype(S))
        new{S,T}()
    end
end

restriction_map(::Type{S}, ::Type{T}) where {S,T} = RestrictionMap{S,T}()
restriction_map(::Type{T}, ::Type{T}) where {T} = IdentityMap{T}()

applymap(map::RestrictionMap{S,T}, x::S) where {S,T} = restrict_space(spacetype(T), x)

left_inverse(map::EmbeddingMap{S,T}) where {S,T} = RestrictionMap{T,S}()
right_inverse(map::RestrictionMap{S,T}) where {S,T} = EmbeddingMap{T,S}()



"A map between isomorphic spaces."
struct IsomorphismMap{S,T} <: AbstractMap{S,T}
    function IsomorphismMap{S,T}() where {S,T}
        @assert embedded(spacetype(T), spacetype(S))
        new{S,T}()
    end
end

isomorphism_map(::Type{S}, ::Type{T}) where {S,T} = IsomorphismMap{S,T}()
isomorphism_map(::Type{T}, ::Type{T}) where {T} = IdentityMap{T}()

applymap(map::IsomorphismMap{S,T}, x::S) where {S,T} = convert_space(spacetype(T), x)

apply_inverse(map::IsomorphismMap{S,T}, y::T) where {S,T} = convert_space(spacetype(S), y)

inv(map::IsomorphismMap{S,T}) where {S,T} = IsomorphismMap{T,S}()
