
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

broadcasted(::DomainSetStyle, ::typeof(+), a::Union{Number,AbstractArray}, d::Domain) =
    map_domain(Translation(a), d)
broadcasted(::DomainSetStyle, ::typeof(+), d::Domain, a::Union{Number,AbstractArray}) =
    map_domain(Translation(a), d)

broadcasted(::DomainSetStyle, ::typeof(-), a::Union{Number,AbstractArray}, d::Domain) =
    map_domain(AffineMap(-1, a), d)
broadcasted(::DomainSetStyle, ::typeof(-), d::Domain, a::Union{Number,AbstractArray}) =
    map_domain(Translation(-a), d)
broadcasted(::DomainSetStyle, ::typeof(-), d::Domain) =
    map_domain(LinearMap(-1), d)

broadcasted(::DomainSetStyle, ::typeof(*), a::Number, d::Domain) =
    map_domain(LinearMap(a), d)
broadcasted(::DomainSetStyle, ::typeof(*), d::Domain, a::Number) =
    map_domain(LinearMap(a), d)


broadcasted(::DomainSetStyle, ::typeof(/), d::Domain, a::Number) =
    mapped_domain(LinearMap(a), d)

broadcasted(::DomainSetStyle, ::typeof(\), a::Number, d::Domain) =
    mapped_domain(LinearMap(a), d)

broadcasted(::DomainSetStyle, m::AbstractMap, d::Domain) = map_domain(m, d)

broadcasted(::DomainSetStyle, fun::Function, d::Domain{T}) where {T} = GenericFunctionMap{T}(fun).(d)
