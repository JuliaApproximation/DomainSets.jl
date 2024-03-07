

"A lazy inverse stores a map `m` and returns `inverse(m, x)`."
struct LazyInverse{T,M} <: SimpleLazyMap{T}
	map	::	M
end

LazyInverse(m) = LazyInverse{codomaintype(m)}(m)
LazyInverse{T}(m) where {T} = LazyInverse{T,typeof(m)}(m)

applymap(m::LazyInverse, x) = inverse(supermap(m), x)

Display.displaystencil(m::LazyInverse) = ["LazyInverse(", supermap(m), ")"]
show(io::IO, mime::MIME"text/plain", m::LazyInverse) = composite_show(io, mime, m)


"""
    inverse(m[, x])

Return the inverse of `m`. The two-argument function evaluates the inverse
at the point `x`.
"""
inverse(m) = LazyInverse(m)
inverse(m::LazyInverse) = supermap(m)
# Concrete maps should implement inverse(m, x)

Base.:\(m::AbstractMap, x) = inverse(m, x)

implements_inverse(m) = !(inverse(m) isa LazyInverse)

"""
    leftinverse(m[, x])

Return a left inverse of the given map. This left inverse `mli` is not unique,
but in any case it is such that `(mli ∘ m) * x = x` for each `x` in the domain
of `m`.

The two-argument function applies the left inverse to the point `x`.
"""
leftinverse(m) = inverse(m)
leftinverse(m, x) = inverse(m, x)

"""
    rightinverse(m[, x])

Return a right inverse of the given map. This right inverse `mri` is not unique,
but in any case it is such that `(m ∘ mri) * y = y` for each `y` in the range
of `m`.

The two-argument function applies the right inverse to the point `x`.
"""
rightinverse(m) = inverse(m)
rightinverse(m, x) = inverse(m, x)

isequalmap(m1::LazyInverse, m2::LazyInverse) = isequalmap(supermap(m1), supermap(m2))
map_hash(m::LazyInverse, h::UInt) = hashrec("LazyInverse", supermap(m), h)
