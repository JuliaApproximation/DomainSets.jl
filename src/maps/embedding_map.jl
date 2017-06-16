# embedding_map.jl

"A map between embedded spaces."
struct EmbeddingMap{T,S} <: AbstractMap{T,S}
    function EmbeddingMap{T,S}() where {T,S}
        @assert embedded(spacetype(S), spacetype(T))
        new{T,S}()
    end
end

embedding_map(::Type{T}, ::Type{S}) where {T,S} = EmbeddingMap{T,S}()
embedding_map(::Type{T}, ::Type{T}) where {T} = IdentityMap{T}()

applymap(map::EmbeddingMap{T,S}, x::S) where {T,S} = convert_space(spacetype(T), x)


"A restriction map from a space to an embedded space."
struct RestrictionMap{T,S} <: AbstractMap{T,S}
    function RestrictionMap{T,S}() where {T,S}
        @assert embedded(spacetype(T), spacetype(S))
        new{T,S}()
    end
end

restriction_map(::Type{T}, ::Type{S}) where {T,S} = RestrictionMap{T,S}()
restriction_map(::Type{T}, ::Type{T}) where {T} = IdentityMap{T}()

applymap(map::RestrictionMap{T,S}, x::S) where {T,S} = restrict_space(spacetype(T), x)



"A map between isomorphic spaces."
struct IsomorphismMap{T,S} <: AbstractMap{T,S}
    function IsomorphismMap{T,S}() where {T,S}
        @assert embedded(spacetype(T), spacetype(S))
        new{T,S}()
    end
end

isomorphism_map(::Type{T}, ::Type{S}) where {T,S} = IsomorphismMap{T,S}()
isomorphism_map(::Type{T}, ::Type{T}) where {T} = IdentityMap{T}()

applymap(map::IsomorphismMap{T,S}, x::S) where {T,S} = convert_space(spacetype(T), x)

apply_inverse(map::IsomorphismMap{T,S}, y::T) where {T,S} = convert_space(spacetype(S), y)

inv(map::IsomorphismMap{T,S}) where {T,S} = IsomorphismMap{S,T}()
