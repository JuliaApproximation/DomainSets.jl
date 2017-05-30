# box.jl

"A BBox is an N-dimensional box specified by its bottom-left and top-right vertices."
struct BBox{N,T}
    left        ::  SVector{N,T}
    right       ::  SVector{N,T}
end

ndims{N,T}(::Type{BBox{N,T}}) = N
ndims{N,T}(::BBox{N,T}) = N

eltype{N,T}(::Type{BBox{N,T}}) = T

# Generic functions for composite types:
element(b::BBox, j::Int) = BBox(b.left[j], b.right[j])
element(b::BBox, range::Range) = BBox(b.left[range], b.right[range])
elements(b::BBox) = tuple([element(b,j) for j in 1:ndims(b)]...)
elements(b::BBox{1}) = (element(b,1),)
elements(b::BBox{2}) = (element(b,1),element(b,2))
elements(b::BBox{3}) = (element(b,1),element(b,2),element(b,3))
elements(b::BBox{4}) = (element(b,1),element(b,2),element(b,3),element(b,4))

nb_elements(b::BBox) = ndims(b)

BBox1{T} = BBox{1,T}
BBox2{T} = BBox{2,T}
BBox3{T} = BBox{3,T}
BBox4{T} = BBox{4,T}

# Dimension-specific constructors
BBox(a::Number, b::Number) = BBox( SVector(a), SVector(b) )
BBox(a, b, c, d) = BBox( SVector(a,c), SVector(b,d) )
BBox(a, b, c, d, e, f) = BBox( SVector(a,c,e), SVector(b,d,f) )
BBox(a, b, c, d, e, f, g, h) = BBox( SVector(a,c,e,g), SVector(b,d,f,h) )
BBox{N,T1,T2}(left::SVector{N,T1},right::SVector{N,T2}) = BBox(SVector{N,promote_type(T1,T2)}(left),SVector{N,promote_type(T1,T2)}(right))

BBox(left, right) = BBox(SVector{length(left)}(left), SVector{length(right)}(right))

convert{T,N}(::Type{BBox{N,T}}, BBox::BBox{N,T}) = BBox
convert{S,T,N}(::Type{BBox{N,S}}, BBox::BBox{N,T}) = BBox{N,S}(SVector{N,S}(left(BBox)),SVector{N,S}(right(BBox)))

tensorproduct(b::BBox) = b
tensorproduct(b::BBox, n::Int) = tensorproduct([b for i=1:n]...)
tensorproduct(b1::BBox, b2::BBox, boxes::BBox...) =
    tensorproduct(BBox(SVector(b1.left..., b2.left...), SVector(b1.right..., b2.right...)), boxes...)

# ⊗{N1,N2,T,S}(b1::BBox{N1,T}, b2::BBox{N2,S}) =
#     BBox(promote_type(T,S)[left(b1)..., left(b2)...], promote_type(T,S)[right(b1)..., right(b2)...])

left(b::BBox) = b.left
left(b::BBox, dim) = b.left[dim]
# In 1D we return a scalar rather than a SVector{1}
left(b::BBox1) = b.left[1]

right(b::BBox) = b.right
right(b::BBox, dim) = b.right[dim]
right(b::BBox1) = b.right[1]

getindex(b::BBox, dim::Int) = (left(b, dim), right(b, dim))

getindex(b::BBox, i::Int, j::Int) = j == 1 ? left(b, i) : right(b, i)

size(b::BBox, dim) = right(b, dim) - left(b, dim)

"Create an equispaced grid on the box with ns[dim] points in each dimension."
equispaced_grid(box, ns) = tensorproduct([PeriodicEquispacedGrid(ns[idx], left(box, idx), right(box, idx)) for idx = 1:ndims(box)]...)

"Create an equispaced grid on the box with n points in each dimension."
equispaced_grid{N}(box::BBox{N}, n::Int) = equispaced_grid(box, ntuple(x->n, Val{N}))

"Create an equispaced grid with n points in the first dimension, and approximately equal stepsize in all dimensions"
equispaced_aspect_grid{N}(box::BBox{N}, n::Int) = equispaced_grid(box, ntuple(i->round(Int,n*(left(box)[i]-right(box)[i])/(left(box)[1]-right(box)[1])), Val{N}))

# Extend a box by a factor of t[i] in each dimension
function extend{N,T}(b::BBox{N,T}, t::SVector{N,T})
    r = SVector{N,T}( [ t[i]*size(b,i) for i in 1:N ] )
    BBox(left(b), left(b) + r)
end

in{N,T}(x, b::BBox{N,T}, dim) = (x[dim] >= left(b,dim)-10eps(T)) && (x[dim] <= right(b,dim)+10eps(T))
in(x, b::BBox) = reduce(&, in(x,b,i) for i = 1:ndims(b))

within(a, b) = (a[1] >= b[1]) && (a[2] <= b[2])
⊂(b1::BBox{1}, b2::BBox{1}) = within(b1[1], b2[1])
⊂(b1::BBox{2}, b2::BBox{2}) = within(b1[1], b2[1]) && within(b1[2], b2[2])
⊂(b1::BBox{3}, b2::BBox{3}) = within(b1[1], b2[1]) && within(b1[2], b2[2]) && within(b1[3], b2[3])
⊂(b1::BBox{4}, b2::BBox{4}) = within(b1[1], b2[1]) && within(b1[2], b2[2]) && within(b1[3], b2[3]) && within(b1[4], b2[4])

## Arithmetic operations

(*)(a::Number,b::BBox) = BBox(a*b.left, a*b.right)
(*)(b::BBox, a::Number) = a*b
(/)(b::BBox, a::Number) = BBox(b.left/a, b.right/a)
(+)(b::BBox, a::AbstractVector) = BBox(b.left+a, b.right+a)
(+)(a::AbstractVector, b::BBox) = b+a
(-)(b::BBox, a::Vector) = BBox(b.left-a, b.right-a)

# Operations on boxes: union and intersection
union(b1::BBox, b2::BBox) = BBox(min.(left(b1),left(b2)), max.(right(b1),right(b2)))

intersect(b1::BBox, b2::BBox) = BBox(max.(left(b1),left(b2)), min.(right(b1),right(b2)))

(+)(b1::BBox, b2::BBox) = union(b1, b2)
(&)(b1::BBox, b2::BBox) = intersect(b1, b2)

volume(box::BBox) = prod(right(box)-left(box))

# There has to be a neater way...
# The implementation of isapprox for Vec is for use in the definition of isapprox for BBox
isapprox(v1::SVector{1}, v2::SVector{1}) = (v1[1] ≈ v2[1])
isapprox(v1::SVector{2}, v2::SVector{2}) = (v1[1] ≈ v2[1]) && (v1[2] ≈ v2[2])
isapprox(v1::SVector{3}, v2::SVector{3}) = (v1[1] ≈ v2[1]) && (v1[2] ≈ v2[2]) && (v1[3] ≈ v2[3])
isapprox(v1::SVector{4}, v2::SVector{4}) = (v1[1] ≈ v2[1]) && (v1[2] ≈ v2[2]) && (v1[3] ≈ v2[3]) && (v1[4] ≈ v2[4])

isapprox(b1::BBox, b2::BBox) = (b1.left ≈ b2.left) && (b1.right ≈ b2.right)

corners(b::BBox{1}) = [b.left[1], b.right[1]]

corners(b::BBox{2}) = [SVector(b.left[1], b.left[2]), SVector(b.left[1], b.right[2]),
    SVector(b.right[1], b.left[2]), SVector(b.right[1], b.right[2])]

function corners{N}(b::BBox{N})
    c1 = corners(element(b,1:N-1))
    c2left = [SVector(c..., b.left[N]) for c in c1]
    c2right = [SVector(c..., b.right[N]) for c in c1]
    vcat(c2left, c2right)
end


show(io::IO, c::BBox{1}) = print(io, "the interval [", left(c, 1), ",", right(c, 1), "]")

show(io::IO, c::BBox{2}) = print(io, "the rectangular box [", left(c, 1), ",", right(c, 1), "] x [", left(c, 2), ",", right(c, 2), "]")

show(io::IO, c::BBox{3}) = print(io, "the cube [", left(c, 1), ",", right(c, 1), "] x [", left(c, 2), ",", right(c, 2), "] x [", left(c, 3), ",", right(c, 3), "]")

show(io::IO, c::BBox{4}) = print(io, "the 4-cube [", left(c, 1), ",", right(c, 1), "] x [", left(c, 2), ",", right(c, 2), "] x [", left(c, 3), ",", right(c, 3), "] x [", left(c, 4), ",", right(c, 4), "]")

# Define the unit box
const unitbox1 = BBox(0, 1)
const unitbox2 = BBox(0, 1, 0, 1)
const unitbox3 = BBox(0, 1, 0, 1, 0, 1)
const unitbox4 = BBox(0, 1, 0, 1, 0, 1, 0, 1)
