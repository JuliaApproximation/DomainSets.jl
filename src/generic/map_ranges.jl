# map_ranges.jl

range(m::AbstractMap) = m * domain(m)

"""
Decide whether a given point `y` is in the range of the map `y=f(x)`, possibly
up to a given tolerance if a tolerance is supplied.
"""
in_range(m::AbstractMap, y) = in(y, range(m))

in_range(m::AbstractMap, y, tolerance) = approx_in(y, range(m), tolerance)


# We can't include domain and ranges of maps in the maps code, because it is
# interpreted before the domain code. Therefore, for a few specific maps, we
# define the domain and range here.

domain(m::IdentityMap{T}) where {T} = FullSpace{T}()
range(m::IdentityMap) = domain(m)

domain(m::ConstantMap{T,S}) where {T,S} = FullSpace{S}()
range(m::ConstantMap) = Point(constant(m))

domain(m::ProductMap) = ProductDomain(map(domain, elements(pm))...)
range(m::ProductMap) = ProductDomain(map(range, elements(pm))...)
