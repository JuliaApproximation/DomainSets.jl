
## Conversion

"""
Convert a domain object to a domain object with the given element type `T`.
"""
convert_domain(::Type{T}, d::Domain{T}) where {T} = d
convert_domain(::Type{T}, d::Domain{S}) where {S,T} = convert(Domain{T}, d)

convert_domain(::Type{T}, d) where {T} = _convert_domain(T, d, eltype(d))
_convert_domain(::Type{T}, d, ::Type{T}) where {T} = d
_convert_domain(::Type{T}, d, ::Type{S}) where {S,T} =
    error("Don't know how to convert the element type of $(typeof(d)) to $S. Please implement `convert_domain` for this combination of types.")

convert_domain(::Type{T}, d::Number) where {T} = convert(T, d)

_convert_domain(::Type{T}, d::AbstractArray{S}, ::Type{S}) where {S,T} = (b = similar(d, T); b .= d)
_convert_domain(::Type{T}, d::Set{S}, ::Type{S}) where {S,T} = convert(Set{T}, d)

## Promotion

"""
Promote the given domains to have a common element type.
"""
promote_domain() = ()
promote_domain(domains...) = promote_domains(domains)

"""
Promote an iterable list of domains to have a common element type.
"""
promote_domains(domains) = _promote_domains(mapreduce(eltype, promote_type, domains), domains)

# We simply rely on Julia's promotion system to promote the element types, and
# we invoke `convert_domain` on the domains with the resulting type.
_promote_domains(::Type{T}, domains) where {T} = convert_domain.(T, domains)
