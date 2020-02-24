
range(m::Map) = m * domain(m)

"""
Decide whether a given point `y` is in the range of the map `y=f(x)`, possibly
up to a given tolerance if a tolerance is supplied.
"""
in_range(m::Map, y) = in(y, range(m))

in_range(m::Map, y, tolerance) = approx_in(y, range(m), tolerance)


# We can't include domain and ranges of maps in the maps code, because it is
# interpreted before the domain code. Therefore, for a few specific maps, we
# define the domain and range here.

domain(m::IdentityMap{T}) where {T} = FullSpace{T}()
range(m::IdentityMap) = domain(m)

domain(m::ConstantMap{T,U}) where {T,U} = FullSpace{T}()
range(m::ConstantMap) = Point(constant(m))

domain(m::ProductMap) = ProductDomain(map(domain, elements(m))...)
range(m::ProductMap) = ProductDomain(map(range, elements(m))...)
