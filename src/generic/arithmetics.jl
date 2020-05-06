# Routines having to do with computations involving domains

+(domain::Domain) = domain

@deprecate *(map::Map, domain::Domain) map.(domain)

promote_numtype(::Type{S}, ::Type{T}) where {S<:Number,T<:Number} = promote_type(S,T)
promote_numtype(::Type{S}, ::Type{SVector{N,T}}) where {N,S<:Number,T<:Number} = SVector{N,promote_type(S,T)}
promote_numtype(::Type{S}, ::Type{T}) where {S<:Number,T} = T

*(a::S, domain::Domain{T}) where {S<:Number,T} =
    map_domain(convert(Map{promote_numtype(S,T)}, a), domain)

*(domain::Domain, a::Number) = a*domain

/(domain::Domain{T}, a::S) where {S<:Number,T} =
    mapped_domain(convert(Map{promote_numtype(S,T)}, a), domain)

@deprecate +(d::Domain, x::Union{AbstractVector,Number}) d .+ x
@deprecate +(x::Union{Number,AbstractVector}, d::Domain) x .+ d
