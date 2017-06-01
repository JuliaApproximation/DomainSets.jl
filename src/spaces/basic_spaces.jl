# basic_spaces.jl

##############################
# Definitions of basic spaces
##############################

"AnySpace is the superset of all possible geometric spaces."
const AnySpace = GeometricSpace{Any}

"A geometric space with integer type."
const IntegerSpace{T <: Integer} = GeometricSpace{T}

"A geometric space with floating point type."
const RealSpace{T <: AbstractFloat} = GeometricSpace{T}

"A geometric space with complex type."
const ComplexPlane{T <: Complex} = GeometricSpace{T}


"The set of integers of type Int (ℤ = \BbbZ)."
const ℤ = IntegerSpace{Int}()

"The set of reals of type Float64 (ℝ = \BbbR)."
const ℝ = RealSpace{Float64}()

"The complex plane with Float64 real and imaginar parts (ℂ = \BbbC)."
const ℂ = ComplexPlane{Complex{Float64}}()

"""
A Euclidean space with static vector types. It is the space of all vectors
of fixed length `N`, with entries of type `T`.
"""
const EuclideanSpace{N,T} = GeometricSpace{SVector{N,T}}

"""
A general Array space, with arrays of dimension `N` and element type `T`.
Note that the arrays can have any size. Thus, `ArraySpace{2,Float64}` contains
all possible matrices of size `m × n` for any combination of `m,n ∈ ℕ`.
"""
const ArraySpace{N,T} = GeometricSpace{Array{T,N}}

^(s::GeometricSpace{T}, ::Val{N}) where {T <: Number,N} = EuclideanSpace{N,T}()

const ℝ2 = ℝ^Val{2}()
const ℝ3 = ℝ^Val{3}()
const ℝ4 = ℝ^Val{4}()

## Isomorphism between scalars and 1D vectors

# Any space with a numeric type can be identified with a 1D euclidean space
embedding_rule(::Type{GeometricSpace{S}}, ::Type{EuclideanSpace{1,T}}) where {S <: Number,T} = _embedding_via_promotion(T,promote_type(S,T))
# and vice-versa
embedding_rule(::Type{EuclideanSpace{1,S}}, ::Type{GeometricSpace{T}}) where {S,T <: Number} = _embedding_via_promotion(T,promote_type(S,T))


# The complex plane can be identified with ℝ^2
embedding_rule(::Type{ComplexPlane{Complex{S}}}, ::Type{EuclideanSpace{2,T}}) where {S <: Number,T} = _embedding_via_promotion(T,promote_type(S,T))
# and vice-versa
embedding_rule(::Type{EuclideanSpace{2,S}}, ::Type{ComplexPlane{Complex{T}}}) where {S, T <: Number} = _embedding_via_promotion(T,promote_type(S,T))
