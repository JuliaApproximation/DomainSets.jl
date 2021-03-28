
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

# Intercept broadcast applied to `in`, e.g. in.(A, d).
# This gives domains an opportunity to provide a more efficient implementation
# when invoked with a set of points `A` at once, especially if `A`
# has particular structure. A common case would be a raster of points,
# for plotting purposes.
# This call can be avoided by typing in.(A, Ref(d)) instead.
broadcasted(::DomainSetStyle, ::typeof(in), A, d::Domain) = broadcast_in(A, d)
broadcasted(::DomainSetStyle, ::typeof(approx_in), A, d::Domain) = broadcast_approx_in(A, d)

"Vectorized version of `in`: apply `x ∈ d` to all elements of `A`."
broadcast_in(A, d::Domain) = in.(A, Ref(d))

"Vectorized version of `approx_in`: apply `x ∈ d` to all elements of `A`."
broadcast_approx_in(A, d::Domain) = approx_in.(A, Ref(d))


## Some arithmetics

# Allow unary minus, but use broadcast for the implementation
-(d::Domain) = (-).(d)

# Allow multiplication and division by numbers, like for vectors
*(a::Number, domain::Domain) = a .* domain
*(domain::Domain, a::Number) = domain .* a
/(domain::Domain, a::Number) = domain ./ a
\(a::Number, domain::Domain) = a .\ domain
