
# We define a new broadcast style for domains, because they may
# represent continuous sets

"The broadcast style associated with domains"
struct DomainSetStyle <: Base.Broadcast.BroadcastStyle end

Base.BroadcastStyle(::Type{<:Domain}) = DomainSetStyle()
Base.BroadcastStyle(::Type{<:AsDomain}) = DomainSetStyle()

Base.broadcastable(d::Domain) = d
Base.broadcastable(d::AsDomain) = d

# DomainSetStyle doesn't mix with ArrayStyle
Base.BroadcastStyle(::DomainSetStyle, ::Base.Broadcast.AbstractArrayStyle) = DomainSetStyle()
Base.BroadcastStyle(::Base.Broadcast.AbstractArrayStyle, ::DomainSetStyle) = DomainSetStyle()

import Base.Broadcast: broadcasted

broadcasted(::DomainSetStyle, ::typeof(+), a::Union{Number,AbstractArray}, d::AnyDomain) =
    map_domain(Translation(a), domain(d))
broadcasted(::DomainSetStyle, ::typeof(+), d::AnyDomain, a::Union{Number,AbstractArray}) =
    map_domain(Translation(a), domain(d))

broadcasted(::DomainSetStyle, ::typeof(-), a::Union{Number,AbstractArray}, d::AnyDomain) =
    map_domain(AffineMap(-1, a), domain(d))
broadcasted(::DomainSetStyle, ::typeof(-), d::AnyDomain, a::Union{Number,AbstractArray}) =
    map_domain(Translation(-a), domain(d))
broadcasted(::DomainSetStyle, ::typeof(-), d::AnyDomain) =
    map_domain(LinearMap{domaineltype(d)}(-1), domain(d))

broadcasted(::DomainSetStyle, ::typeof(*), a::Number, d::AnyDomain) =
    map_domain(LinearMap{domaineltype(d)}(a), domain(d))
broadcasted(::DomainSetStyle, ::typeof(*), d::AnyDomain, a::Number) =
    map_domain(LinearMap{domaineltype(d)}(a), domain(d))


broadcasted(::DomainSetStyle, ::typeof(/), d::AnyDomain, a::Number) =
    mapped_domain(LinearMap(a), domain(d))

broadcasted(::DomainSetStyle, ::typeof(\), a::Number, d::AnyDomain) =
    mapped_domain(LinearMap(a), domain(d))

broadcasted(::DomainSetStyle, m::AbstractMap, d::AnyDomain) =
    map_domain(m, domain(d))

broadcasted(::DomainSetStyle, fun::Function, d::AnyDomain) =
    convert(Map{domaineltype(d)}, fun).(d)

# Intercept broadcast applied to `in`, e.g. in.(A, d).
# This gives domains an opportunity to provide a more efficient implementation
# when invoked with a set of points `A` at once, especially if `A`
# has particular structure. A common case would be a raster of points,
# for plotting purposes.
# This call can be avoided by typing in.(A, Ref(d)) instead.
broadcasted(::DomainSetStyle, ::typeof(in), A, d::AnyDomain) =
    vectorized_in(A, domain(d))
broadcasted(::DomainSetStyle, ::typeof(approx_in), A, d::AnyDomain) =
    vectorized_approx_in(A, domain(d))
broadcasted(::DomainSetStyle, ::typeof(approx_in), A, d::AnyDomain, tol) =
    vectorized_approx_in(A, domain(d), tol)

broadcasted(::DomainSetStyle, ::typeof(∉), A, d::AnyDomain) = A .∉ Ref(domain(d))

@deprecate broadcast_in(A, d::Domain) vectorized_in(A, d)
@deprecate broadcast_approx_in(A, d::Domain) vectorized_approx_in(A, d)
@deprecate broadcast_approx_in(A, d::Domain, tol) vectorized_approx_in(A, d, tol)

"Vectorized version of `in`: apply `x ∈ d` to all elements of `A`."
vectorized_in(A, d) = in.(A, Ref(d))

"Vectorized version of `approx_in`: apply `x ∈ d` to all elements of `A`."
vectorized_approx_in(A, d) = approx_in.(A, Ref(d))
vectorized_approx_in(A, d, tol) = approx_in.(A, Ref(d), tol)


## Some arithmetics

# Allow unary minus, but use broadcast for the implementation
-(d::AnyDomain) = (-).(d)

# Allow multiplication and division by numbers, like for vectors
*(a::Number, domain::AnyDomain) = a .* domain
*(domain::AnyDomain, a::Number) = domain .* a
/(domain::AnyDomain, a::Number) = domain ./ a
\(a::Number, domain::AnyDomain) = a .\ domain
