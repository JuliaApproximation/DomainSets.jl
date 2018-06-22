# basic_spaces.jl

##############################
# Definitions of basic spaces
##############################

"A geometric space with integer type."
const IntegerSpace{T <: Integer} = GeometricSpace{T}

"A geometric space with rational types."
const RationalSpace{T} = GeometricSpace{Rational{T}}

subeltype(::Type{RationalSpace{T}}) where {T} = T
similar_space(::Type{RationalSpace{T}}, ::Type{S}) where {T,S} = RationalSpace{S}

"A geometric space with floating point type."
const RealSpace{T <: AbstractFloat} = GeometricSpace{T}

"A geometric space with complex type."
const ComplexSpace{T} = GeometricSpace{Complex{T}}

subeltype(::Type{ComplexSpace{T}}) where {T} = T
similar_space(::Type{ComplexSpace{T}}, ::Type{S}) where {T,S} = ComplexSpace{S}


"""
A Euclidean space with static vector types. It is the space of all vectors
of fixed length `N`, with entries of type `T`.
"""
const VectorSpace{N,T} = GeometricSpace{SVector{N,T}}
const EuclideanSpace{N,T} = VectorSpace{N,T}

subeltype(::Type{VectorSpace{N,T}}) where {N,T} = T
similar_space(::Type{VectorSpace{N,T}}, ::Type{S}) where {N,T,S} = VectorSpace{N,S}

"""
A general Array space, with arrays of dimension `N` and element type `T`.
Note that the arrays can have any size. Thus, `ArraySpace{2,Float64}` contains
all possible matrices of size `m × n` for any combination of `m,n ∈ ℕ`.
"""
const ArraySpace{N,T} = GeometricSpace{Array{T,N}}

subeltype(::Type{ArraySpace{N,T}}) where {N,T} = T
similar_space(::Type{ArraySpace{N,T}}, ::Type{S}) where {N,T,S} = ArraySpace{N,S}

## Convenience functions

^(::Type{GeometricSpace{T}}, ::Type{Val{N}}) where {T <: Number,N} = VectorSpace{N,T}

## Some standard spaces

"The set of integers of type Int (ℤ = \\BbbZ)."
const ℤ = IntegerSpace{Int}
"The set of rational numbers of type Rational{Int} (ℚ = \\BbbQ)."
const ℚ = RationalSpace{Int}
"The set of reals of type Float64 (ℝ = \\BbbR)."
const ℝ = RealSpace{Float64}
"The complex plane with Float64 real and imaginar parts (ℂ = \\BbbC)."
const ℂ = ComplexSpace{Float64}

"The space ℝ^1"
const ℝ1 = ℝ^Val{1}
"The space ℝ^2"
const ℝ2 = ℝ^Val{2}
"The space ℝ^3"
const ℝ3 = ℝ^Val{3}
"The space ℝ^4"
const ℝ4 = ℝ^Val{4}




##############################
# Embeddings and conversions
##############################

## Isomorphism between scalars and 1D vectors

isomorphism_reduction(::Type{GeometricSpace{T}}, ::Type{VectorSpace{1,S}}) where {T <: Number,S} =
    (GeometricSpace{T},GeometricSpace{S})

convert_space(::Type{GeometricSpace{T}}, x::SVector{1,T}) where {T} = x[1]
convert_space(::Type{VectorSpace{1,T}}, x::T) where {T} = SVector(x)

## Isomorphism between the complex plane and ℝ^2

isomorphism_reduction(::Type{VectorSpace{2,T}}, ::Type{ComplexSpace{S}}) where {T,S} =
    (GeometricSpace{T},GeometricSpace{S})

convert_space(::Type{ComplexSpace{T}}, x::SVector{2,T}) where {T} = x[1] + im*x[2]
convert_space(::Type{VectorSpace{2,T}}, x::Complex{T}) where {T} = SVector(real(x), imag(x))

## Isomorphism between vector spaces of the same dimension

# Normally, this line never gets called when T equals S, because that case is
# caught higher up in the algorithms.
isomorphism_reduction(::Type{VectorSpace{N,T}}, ::Type{VectorSpace{N,S}}) where {N,T,S} =
    (GeometricSpace{T},GeometricSpace{S})

# Just in case, we make sure if it ever happens that we notice
isomorphism_reduction(::Type{VectorSpace{N,T}}, ::Type{VectorSpace{N,T}}) where {N,T} =
    error("isomorphism_reduction was called with identical types")

convert_space(A::Type{VectorSpace{N,T}}, x::SVector{N,T}) where {N,T} = x
convert_space(A::Type{VectorSpace{N,T}}, x::SVector{N,S}) where {N,T,S} =
    SVector{N,T}(convert_space.(GeometricSpace{T}, x))

restrict_space(A::Type{VectorSpace{N,T}}, x::SVector{N,T}) where {N,T} = x
restrict_space(A::Type{VectorSpace{N,T}}, x::SVector{N,S}) where {N,T,S} =
    SVector{N,T}(restrict_space.(GeometricSpace{T}, x))


## For completeness, we have an isomorphism rule for ComplexSpace and RationalSpace
## to themselves as well. These are probably not useful.

isomorphism_reduction(::Type{ComplexSpace{T}}, ::Type{ComplexSpace{S}}) where {T,S} =
    (GeometricSpace{T},GeometricSpace{S})

convert_space(A::Type{ComplexSpace{T}}, x::Complex{T}) where {T} = x
convert_space(A::Type{ComplexSpace{T}}, x::Complex{S}) where {T,S} =
    convert_space(GeometricSpace{T}, real(x)) + im * convert_space(GeometricSpace{T}, imag(x))

restrict_space(A::Type{ComplexSpace{T}}, x::Complex{T}) where {T} = x
restrict_space(A::Type{ComplexSpace{T}}, x::Complex{S}) where {T,S} =
    restrict_space(GeometricSpace{T}, real(x)) + im * restrict_space(GeometricSpace{T}, imag(x))

isomorphism_reduction(::Type{RationalSpace{T}}, ::Type{RationalSpace{S}}) where {T,S} =
    (GeometricSpace{T},GeometricSpace{S})

convert_space(A::Type{RationalSpace{T}}, x::Rational{T}) where {T} = x
convert_space(A::Type{RationalSpace{T}}, x::Rational{S}) where {T,S} =
    convert_space(GeometricSpace{T}, numerator(x)) // convert_space(GeometricSpace{T}, denominator(x))

restrict_space(A::Type{RationalSpace{T}}, x::Rational{T}) where {T} = x
restrict_space(A::Type{RationalSpace{T}}, x::Rational{S}) where {T,S} =
    restrict_space(GeometricSpace{T}, numerator(x)) // restrict_space(GeometricSpace{T}, denominator(x))

## Embedding of R^N in R^(N+1)

# We only go up to dimension four
embedding_reduction(::Type{VectorSpace{1,T}}, ::Type{VectorSpace{2,S}}) where {T,S} =
    (GeometricSpace{T},GeometricSpace{S})
embedding_reduction(::Type{VectorSpace{2,T}}, ::Type{VectorSpace{3,S}}) where {T,S} =
    (GeometricSpace{T},GeometricSpace{S})
embedding_reduction(::Type{VectorSpace{3,T}}, ::Type{VectorSpace{4,S}}) where {T,S} =
    (GeometricSpace{T},GeometricSpace{S})

convert_space(::Type{VectorSpace{2,T}}, x::SVector{1,T}) where {T} = SVector(x[1], 0)
convert_space(::Type{VectorSpace{3,T}}, x::SVector{2,T}) where {T} = SVector(x[1], x[2], 0)
convert_space(::Type{VectorSpace{4,T}}, x::SVector{3,T}) where {T} = SVector(x[1], x[2], x[3], 0)

# Since these are true embeddings, we also provide a left inverse
restrict_space(::Type{VectorSpace{1,T}}, x::SVector{2,T}) where {T} = SVector(x[1])
restrict_space(::Type{VectorSpace{2,T}}, x::SVector{3,T}) where {T} = SVector(x[1], x[2])
restrict_space(::Type{VectorSpace{3,T}}, x::SVector{4,T}) where {T} = SVector(x[1], x[2], x[3])

# Superspaces of a vector space have one higher dimension
superspace(::Type{VectorSpace{1,T}}) where {T} = VectorSpace{2,T}
superspace(::Type{VectorSpace{2,T}}) where {T} = VectorSpace{3,T}
superspace(::Type{VectorSpace{3,T}}) where {T} = VectorSpace{4,T}
