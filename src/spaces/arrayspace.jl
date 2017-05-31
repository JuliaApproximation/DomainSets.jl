# arrayspace.jl

"""
An `ArraySpace` is a space where each element is an array of dimension `N` and
with eltype `T`. The size of the array is not part of the type and is stored
in a field of a concrete instance.
"""
struct ArraySpace{T,N} <: GeometricSpace{Array{T,N}}
    size    ::  NTuple{N,Int}
end

const VectorSpace{T} = ArraySpace{T,1}
const MatrixSpace{T} = ArraySpace{T,2}
const TensorSpace{T} = ArraySpace{T,3}

# This constructor deduces N from the given tuple
ArraySpace{T}(s::NTuple{N,Int}) where {T,N} = ArraySpace{T,N}(s)

# We convert a list of Int's to a tuple
# - with N not specified
ArraySpace{T}(s::Int...) where {T} = ArraySpace{T}(s)
# - with N specified
ArraySpace{T,N}(s::Int...) where {T,N} = ArraySpace{T,N}(s)

# Finally, we use T = Float64 as default type if not specified
ArraySpace(args...) = ArraySpace{Float64}(args...)
VectorSpace(args...) = VectorSpace{Float64}(args...)
MatrixSpace(args...) = MatrixSpace{Float64}(args...)
TensorSpace(args...) = TensorSpace{Float64}(args...)

size(space::ArraySpace) = space.size

array_eltype{T,N}(space::ArraySpace{T,N}) = T

zero{T,N}(space::ArraySpace{T,N}) = zeros(T, size(space))

in{T,N}(x::Array{T,N}, space::ArraySpace{T,N}) = size(x) == size(space)
