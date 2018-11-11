
## Conversion

convert_domain(::Type{T}, d::Domain{T}) where {T} = d
convert_domain(::Type{S}, d::Domain{T}) where {S,T} = convert(Domain{S}, d)

convert_domain(::Type{T}, d) where {T} = _convert_domain(T, d, eltype(d))
_convert_domain(::Type{T}, d, ::Type{T}) where {T} = d
_convert_domain(::Type{T}, d, ::Type{S}) where {S,T} =
    error("Don't know how to convert the element type of $(typeof(d)) to $S. Please implement `convert_domain` for this combination of types.")

convert_domain(::Type{T}, d::Number) where {T} = convert(T, d)

convert_domain(::Type{T}, d::AbstractArray{T}) where {T} = d
convert_domain(::Type{T}, d::AbstractArray{S}) where {S,T} = (b = similar(d, T); b .= d)


## Promotion

# We want to make sure that a set of domains, that may have different types,
# are promoted to domains with a joined element type.
# We simply rely on Julia's promotion system to promote the element types, and
# we invoke `convert_domain` on the domains with the resulting type.

promote_domain() = ()
promote_domain(a) = (a,)

promote_domain(x, y) = _promote_domain(promote_type(eltype(x), eltype(y)), x, y)
promote_domain(x, y, zs...) = _promote_domain(promote_type(map(eltype, (x, y, zs...))), x, y, zs...)

_promote_domain(::Type{T}, x, y) where {T} = (convert_domain(T, x), convert_domain(T, y))
_promote_domain(::Type{T}, x, y, zs...) where {T} =
    map(d -> convert_domain(T, d), (x, y, zs...))
