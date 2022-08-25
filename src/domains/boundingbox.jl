
boundingbox(d::Vector{T}) where {T <: Number} = minimum(d)..maximum(d)
boundingbox(d::Set{T}) where {T<:Number} = minimum(d)..maximum(d)

"Return the bounding box of the union of two or more bounding boxes."
unionbox(d::Domain) = d
unionbox(d1::Domain, d2::Domain) = unionbox(promote_domains(d1, d2)...)
unionbox(d1::Domain, d2::Domain, domains...) =
    unionbox(unionbox(d1,d2), domains...)

unionbox(d1::Domain{T}, d2::Domain{T}) where {T} = unionbox1(d1, d2)
unionbox1(d1, d2) = unionbox2(d1, d2)
unionbox2(d1, d2) = FullSpace{eltype(d1)}()
unionbox1(d1::EmptySpace, d2) = d2
unionbox1(d1::FullSpace, d2) = d1
unionbox2(d1, d2::EmptySpace) = d1
unionbox2(d1, d2::FullSpace) = d2

unionbox(d1::D, d2::D) where {D<:FixedInterval} = d1

function unionbox(d1::AbstractInterval{T}, d2::AbstractInterval{T}) where {T}
    a, b = endpoints(d1)
    c, d = endpoints(d2)
    A = min(a, c)
    B = max(b, d)
    isinf(A) && isinf(B) ? FullSpace{T}() : A..B
end

unionbox(d1::HyperRectangle{T}, d2::HyperRectangle{T}) where {T} =
    Rectangle{T}(map(unionbox, components(d1), components(d2)))

"Return the bounding box of the intersection of two or more bounding boxes."
intersectbox(d::Domain) = d
intersectbox(d1::Domain, d2::Domain) = intersectbox(promote_domains(d1, d2)...)
intersectbox(d1::Domain, d2::Domain, domains...) =
    intersectbox(intersectbox(d1,d2), domains...)

intersectbox(d1::Domain{T}, d2::Domain{T}) where {T} = intersectbox1(d1, d2)
intersectbox1(d1, d2) = intersectbox2(d1, d2)
intersectbox2(d1, d2) = FullSpace{eltype(d1)}()
intersectbox1(d1::EmptySpace, d2) = d1
intersectbox1(d1::FullSpace, d2) = d2
intersectbox2(d1, d2::EmptySpace) = d2
intersectbox2(d1, d2::FullSpace) = d1

intersectbox(d1::D, d2::D) where {D<:FixedInterval} = d1

intersectbox(d1::AbstractInterval{T}, d2::AbstractInterval{T}) where {T} =
    intersectdomain(d1, d2)

function intersectbox(d1::HyperRectangle{T}, d2::HyperRectangle{T}) where {T}
    d = Rectangle{T}(map(intersectbox, components(d1), components(d2)))
    isempty(d) ? EmptySpace{T}() : d
end

boundingbox(d::AbstractMappedDomain) = map_boundingbox(boundingbox(superdomain(d)), forward_map(d))

function map_boundingbox(box::AbstractInterval, fmap)
    l,r = (leftendpoint(box),rightendpoint(box))
    ml = fmap(l); mr = fmap(r)
    min(ml,mr)..max(ml,mr)
end

# This is a best effort implementation, it could be wrong for some maps,
# because we only map the corners. Hence we restrict to affine maps here.
map_boundingbox(box::HyperRectangle, fmap::AbstractAffineMap) =
    map_boundingbox_generic(box, fmap)

# This is a best effort implementation. It could be wrong for some maps,
# because we only map the corners.
function map_boundingbox_generic(box::HyperRectangle{T}, fmap) where {T}
    mapped_corners = map(fmap, corners(box))
    left = [minimum(x[j] for x in mapped_corners) for j in 1:dimension(box)]
    right = [maximum(x[j] for x in mapped_corners) for j in 1:dimension(box)]
    Rectangle{T}(left, right)
end
