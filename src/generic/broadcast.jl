
# We define a new broadcast style for domains, because they may
# represent continuous sets

"The broadcast style associated with domains"
struct DomainSetStyle <: Base.Broadcast.BroadcastStyle end

Base.BroadcastStyle(::Type{<:Domain}) = DomainSetStyle()

Base.broadcastable(d::Domain) = d

# DomainSetStyle doesn't mix with ArrayStyle
Base.BroadcastStyle(::DomainSetStyle, ::Base.Broadcast.AbstractArrayStyle) = DomainSetStyle()
Base.BroadcastStyle(::Base.Broadcast.AbstractArrayStyle, ::DomainSetStyle) = DomainSetStyle()

import Base.Broadcast: broadcasted

broadcasted(::DomainSetStyle, ::typeof(+), a::Number, d::Domain{T}) where {T} =
    map_domain(Translation{promote_type(T,typeof(a))}(a), d)
broadcasted(::DomainSetStyle, ::typeof(+), d::Domain{T}, a::Number) where {T} =
    map_domain(Translation{promote_type(T,typeof(a))}(a), d)

broadcasted(::DomainSetStyle, ::typeof(+), a::AbstractArray{S}, d::EuclideanDomain{N,T}) where {N,S,T} =
    map_domain(Translation(SVector{N,promote_type(S,T)}(a)), d)
broadcasted(::DomainSetStyle, ::typeof(+), d::EuclideanDomain{N,T}, a::AbstractArray{S}) where {N,S,T} =
    map_domain(Translation(SVector{N,promote_type(S,T)}(a)), d)

broadcasted(::DomainSetStyle, m::AbstractMap, d::Domain) = map_domain(m, d)
