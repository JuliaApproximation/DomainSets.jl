# euclidean.jl

"""
A euclidean space is a space of the form T^N, where T is an element type and
N is the dimension.
"""
struct EuclideanSpace{N,T} <: GeometricSpace{SVector{N,T}}
end

const ESpace1d{T} = EuclideanSpace{1,T}
const ESpace2d{T} = EuclideanSpace{2,T}
const ESpace3d{T} = EuclideanSpace{3,T}
const ESpace4d{T} = EuclideanSpace{4,T}

ndims{N,T}(::EuclideanSpace{N,T}) = N



## Some arithmetic with spaces

# These are not type-stable, but convenient:
^{T}(space::UnivariateSpace{T}, n::Int) = EuclideanSpace{n,T}()
^{T,N}(space::EuclideanSpace{N,T}, n::Int) = EuclideanSpace{N*n,T}()


## Isomorphisms

# A real 2d Euclidean space is isomorphic to the complex plane
isomorphic{T <: Real}(space1::EuclideanSpace{2,T}, space2::ComplexPlane{T}) = true
isomorphic{T <: Real}(space1::ComplexPlane{T}, space2::EuclideanSpace{2,T}) = true

# Convert back and forth using real and imaginary parts
unsafe_convert_space(x, space1::EuclideanSpace{2}, space2::ComplexPlane) = x[1] + im*x[2]
unsafe_convert_space(x, space1::ComplexPlane, space2::EuclideanSpace{2}) = SVector(real(x), imag(x))

## Embedding declarations

# Euclidean spaces with the same type T are embedded if N1 ≦ N2
embedded{N1,N2,T}(s1::EuclideanSpace{N1,T}, s2::EuclideanSpace{N2,T}) = N1 <= N2

# Euclidean spaces with the same dimension are embedded if T1 can be promoted to T2
embedded{N,T1,T2}(s1::EuclideanSpace{N,T1}, s2::EuclideanSpace{N,T2}) = promote_type(T1,T2) == T2


## Embedding logic

# This is a general implementation, that is probably slow
unsafe_promote_space{N1,N2,T}(x, s1::EuclideanSpace{N1,T}, s2::EuclideanSpace{N2,T}) =
    SVector{N2,T}(x..., zeros(SVector{T,N2-N1})...)

# Hence, we provide some shortcuts
# - Extend a point from 1 to 2 dimensions
unsafe_promote_space{T}(x, s1::ESpace1d{T}, s2::ESpace2d{T}) =
    SVector{2,T}(x[1], 0)

# - Extend a point from 1 to 3 dimensions
unsafe_promote_space{T}(x, s1::ESpace1d{T}, s2::ESpace3d{T}) =
    SVector{3,T}(x[1], 0, 0)

# - Extend a point from 2 to 3 dimensions
unsafe_promote_space{T}(x, s1::ESpace2d{T}, s2::ESpace3d{T}) =
    SVector{3,T}(x[1], x[2], 0)

# - Promote the numeric type of a point
unsafe_promote_space{N,T1,T2}(x, s1::EuclideanSpace{N,T1}, s2::EuclideanSpace{N,T2}) =
    SVector{N,T2}(x)
