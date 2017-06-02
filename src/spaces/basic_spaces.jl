# basic_spaces.jl

##############################
# Definitions of basic spaces
##############################

"A geometric space with integer type."
const IntegerSpace{T <: Integer} = GeometricSpace{T}

"A geometric space with rational types."
const RationalSpace{T <: Rational} = GeometricSpace{T}

"A geometric space with floating point type."
const RealSpace{T <: AbstractFloat} = GeometricSpace{T}

"A geometric space with complex type."
const ComplexPlane{T <: Complex} = GeometricSpace{T}

"""
A Euclidean space with static vector types. It is the space of all vectors
of fixed length `N`, with entries of type `T`.
"""
const VectorSpace{N,T} = GeometricSpace{SVector{N,T}}
const EuclideanSpace{N,T} = VectorSpace{N,T}

"""
A general Array space, with arrays of dimension `N` and element type `T`.
Note that the arrays can have any size. Thus, `ArraySpace{2,Float64}` contains
all possible matrices of size `m × n` for any combination of `m,n ∈ ℕ`.
"""
const ArraySpace{N,T} = GeometricSpace{Array{T,N}}

## Convenience functions

^(::GeometricSpace{T}, ::Val{N}) where {T <: Number,N} = VectorSpace{N,T}()

## Some standard spaces

"The set of integers of type Int (ℤ = \BbbZ)."
const ℤ = IntegerSpace{Int}()
"The set of rational numbers of type Rational{Int} (ℚ = \BbbQ)."
const ℚ = RationalSpace{Rational{Int}}()
"The set of reals of type Float64 (ℝ = \BbbR)."
const ℝ = RealSpace{Float64}()
"The complex plane with Float64 real and imaginar parts (ℂ = \BbbC)."
const ℂ = ComplexPlane{Complex{Float64}}()

"The space ℝ^2"
const ℝ2 = ℝ^Val{2}()
"The space ℝ^3"
const ℝ3 = ℝ^Val{3}()
"The space ℝ^4"
const ℝ4 = ℝ^Val{4}()




##############################
# Embeddings and conversions
##############################

## Isomorphism between scalars and 1D vectors

# Any space with a numeric type can be identified with a 1D euclidean space
embedding(::Type{GeometricSpace{S}}, ::Type{VectorSpace{1,T}}) where {S <: Number,T} = promotes_to(S,T)
# and vice-versa
embedding(::Type{VectorSpace{1,S}}, ::Type{GeometricSpace{T}}) where {S,T <: Number} = promotes_to(S,T)

convert_space(::Type{GeometricSpace{T}}, x::SVector{1,S}) where {T <: Number,S} = T(x[1])
convert_space(::Type{VectorSpace{1,T}}, x::S) where {T,S <: Number} = SVector{1,T}(x)

## Isomorphism between the complex plane and ℝ^2

# The complex plane is embedded in ℝ^2
embedding(::Type{ComplexPlane{Complex{S}}}, ::Type{VectorSpace{2,T}}) where {S,T} = promotes_to(S,T)
# and vice-versa
embedding(::Type{VectorSpace{2,S}}, ::Type{ComplexPlane{Complex{T}}}) where {S,T} = promotes_to(S,T)

convert_space(::Type{ComplexPlane{Complex{T}}}, x::SVector{2,S}) where {T,S} = T(x[1]) + im*T(x[2])
convert_space(::Type{VectorSpace{2,T}}, x::Complex{S}) where {T,S} = SVector{2,T}(real(x), imag(x))

# A smaller Euclidean space can be embedded into a larger one. The canonical
# embedding we implement is extension by zero.
# In general, this calls for generated functions. We include some specific rules
# here in low dimensions (up to four).
embedding(::Type{VectorSpace{1,S}}, ::Type{VectorSpace{2,T}}) where {S,T} = promotes_to(S,T)
embedding(::Type{VectorSpace{1,S}}, ::Type{VectorSpace{4,T}}) where {S,T} = promotes_to(S,T)
embedding(::Type{VectorSpace{2,S}}, ::Type{VectorSpace{3,T}}) where {S,T} = promotes_to(S,T)
embedding(::Type{VectorSpace{3,S}}, ::Type{VectorSpace{4,T}}) where {S,T} = promotes_to(S,T)

convert_space(::Type{VectorSpace{2,T}}, x::SVector{1,S}) where {T,S} = SVector{2,T}(x[1], 0)
convert_space(::Type{VectorSpace{3,T}}, x::SVector{2,S}) where {T,S} = SVector{3,T}(x[1], x[2], 0)
convert_space(::Type{VectorSpace{4,T}}, x::SVector{3,S}) where {T,S} = SVector{4,T}(x[1], x[2], x[3], 0)


##############
# Promotions
##############

# promote any numeric type to floating point
superspace(::Type{GeometricSpace{T}}) where {T <: Number} = RealSpace{float(T)}
# promote floating point to a 1D-vector
superspace(::Type{GeometricSpace{T}}) where {T <: AbstractFloat} = VectorSpace{1,T}

# promote the complex plane to ℝ^2
superspace(::Type{ComplexPlane{Complex{T}}}) where {T} = VectorSpace{2,T}

# promote a vector space to a space of one more dimension
superspace(::Type{VectorSpace{1,T}}) where {T} = VectorSpace{2,T}
superspace(::Type{VectorSpace{2,T}}) where {T} = VectorSpace{3,T}
superspace(::Type{VectorSpace{3,T}}) where {T} = VectorSpace{4,T}

# favour ℝ^2 over ℂ, if the type T is the same (otherwise the space with the widest T is chosen by the default rules)
isomorphism_promotion_rule(::Type{ComplexPlane{Complex{T}}}, ::Type{VectorSpace{2,T}}) where {T} = VectorSpace{2,T}
