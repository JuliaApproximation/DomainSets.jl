# mapped_domain.jl

##################
# Mapped domains
##################

"""
A MappedDomain consists of a domain and a bidirectional map. The forward map
maps the domain onto the mapped domain, and the inverse map maps it back.
A point lies in the mapped domain, if the inverse map of that point lies in the
original domain.
"""
# TODO: experiment with leaving out the type parameters and implement fast indomain_grid
struct MappedDomain{DOMAIN <: Domain,MAP,T} <: Domain{T}
    domain  ::  DOMAIN
    # The forward map, from the underlying domain to the mapped domain
    fmap    ::  MAP

    # With this inner constructor we enforce that N is the dimension of the domain
    MappedDomain{DOMAIN,MAP,T}(domain::Domain{T}, fmap) where {DOMAIN <: Domain,MAP,T} = new{DOMAIN,MAP,T}(domain, fmap)
end

MappedDomain(domain::Domain{T}, fmap) where {T} =
    MappedDomain{typeof(domain),typeof(fmap),T}(domain, fmap)

domain(d::MappedDomain) = d.domain

mapping(d::MappedDomain) = d.fmap

indomain(x, d::MappedDomain) = indomain(inverse_map(mapping(d), x), domain(d))

# Now here is a problem: how do we compute a bounding box, without extra knowledge
# of the map? We can only do this for some maps.
boundingbox(d::MappedDomain) = mapped_boundingbox(boundingbox(domain(d)), mapping(d))

function mapped_boundingbox(box::BBox1, fmap)
    l,r = box[1]
    ml = fmap*l
    mr = fmap*r
    BBox(min(ml,mr), max(ml,mr))
end

# In general, we can at least map all the corners of the bounding box of the
# underlying domain, and compute a bounding box for those points. This will be
# correct for affine maps.
function mapped_boundingbox{N}(box::BBox{N}, fmap)
    crn = corners(box)
    mapped_corners = [fmap*c for c in crn]
    left = [minimum([mapped_corners[i][j] for i in 1:length(mapped_corners)]) for j in 1:N]
    right = [maximum([mapped_corners[i][j] for i in 1:length(mapped_corners)]) for j in 1:N]
    BBox(left, right)
end

# We can do better for diagonal maps, since the problem simplifies: each dimension
# is mapped independently.
mapped_boundingbox{N}(box::BBox{N}, fmap::ProductMap) =
    tensorproduct([mapped_boundingbox(element(box,i), element(fmap,i)) for i in 1:N]...)

apply_map(domain::Domain, map::AbstractMap) = MappedDomain(domain, map)

apply_map(d::MappedDomain, map::AbstractMap) = MappedDomain(domain(d), map*mapping(d))

(*)(map::AbstractMap, domain::Domain) = apply_map(domain, map)

(*)(domain::Domain, a::Number) = scaling_map(a*diagm(ones(eltype(domain)))) * domain

# TODO: revise
(+){N,T}(d::Domain{N}, x::SVector{N,T}) = AffineMap(eye(SMatrix{N,N,T}),x) * d
(+){N}(d::Domain{N}, x::AbstractVector) = d + SVector{N}(x)

show(io::IO, d::MappedDomain) =  print(io, "A mapped domain based on ", domain(d))
