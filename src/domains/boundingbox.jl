
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
    Rectangle{T}(map(unionbox, elements(d1), elements(d2)))

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
    d = Rectangle{T}(map(intersectbox, elements(d1), elements(d2)))
    isempty(d) ? EmptySpace{T}() : d
end
