# arithmetics.jl
# Routines having to do with computations involving domains

(*)(map::AbstractMap, domain::Domain) = applymap(map, domain)

(*)(domain::Domain, a::Number) = scaling_map(a*diagm(ones(eltype(domain)))) * domain

# TODO: revise
(+)(d::Domain, x::SVector{N,T}) where {N,T} = AffineMap(eye(SMatrix{N,N,T}),x) * d
# (+){N}(d::Domain{N}, x::AbstractVector) = d + SVector{N}(x)
