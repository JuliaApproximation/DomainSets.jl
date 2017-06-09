# arithmetics.jl
# Routines having to do with computations involving domains

(*)(map::AbstractMap, domain::Domain) = map_domain(map, domain)

(*)(a::Number, domain::Domain{T}) where {T} = LinearMap{T}(1/a) * domain
(*)(domain::Domain, a::Number) = a*domain

(/)(domain::Domain{T}, a::Number) where {T} = LinearMap{T}(a) * domain

(+)(d::Domain, x::SVector{N,T}) where {N,T} = Translation(-x) * d
(+)(x::SVector, d::Domain) = d+x
