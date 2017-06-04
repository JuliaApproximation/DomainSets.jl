# productspace.jl

"""
By convention a geometric space with a tuple type `T` represents a product space
of the individual entries of `T`.
"""
const ProductSpace{T <: Tuple} = GeometricSpace{T}

tensorproduct(a::GeometricSpace) = tensorproduct(typeof(a))
tensorproduct(a::GeometricSpace, b::GeometricSpace) = tensorproduct(typeof(a), typeof(b))()
tensorproduct(a::GeometricSpace, b::GeometricSpace, c::GeometricSpace) = tensorproduct(typeof(a), typeof(b), typeof(c))()
tensorproduct(a::GeometricSpace, b::GeometricSpace, c::GeometricSpace, d::GeometricSpace) = tensorproduct(typeof(a), typeof(b), typeof(c), typeof(d))()

tensorproduct(::Type{GeometricSpace{T}}) where {T} = ProductSpace{Tuple{T}}
tensorproduct(::Type{GeometricSpace{T}}, ::Type{GeometricSpace{S}}) where {T,S} = ProductSpace{Tuple{T,S}}
tensorproduct(::Type{GeometricSpace{T}}, ::Type{GeometricSpace{S}}, ::Type{GeometricSpace{U}}) where {T,S,U} = ProductSpace{Tuple{T,S,U}}
tensorproduct(::Type{GeometricSpace{T}}, ::Type{GeometricSpace{S}}, ::Type{GeometricSpace{U}}, ::Type{GeometricSpace{V}}) where {T,S,U,V} = ProductSpace{Tuple{T,S,U,V}}

## Superspaces

# Tensor product spaces may be nested, and we intend to flatten them. We proceed
# by removing tuples one at a time. We do this manually, and only for cases up
# to four dimensions in total.
#
# The most natural isomorphism is between tuples with homogeneous types, and
# static vectors: (T,T,T) and [T,T,T].
# A nested example is (T,(T,T)), which we identify with (T,T,T).
# A more nested case is (T,(T,(T,))), which we identify with (T,(T,T)), and so on.
# We try to be general but focus on numeric types T in most cases,
# i.e. we restrict to T <: Number, because T might otherwise match a Tuple.

# The correspondence between (T,T,...T) and [T,T,T] is a one-liner.
# Up to 4D this covers (T,), (T,T), (T,T,T), (T,T,T,T)
superspace(::Type{GeometricSpace{NTuple{N,T}}}) where {N,T <: Number} = VectorSpace{N,T}

# In 1D we can have nested tuples ((T,),)
superspace(::Type{GeometricSpace{Tuple{Tuple{T}}}}) where {T} = GeometricSpace{Tuple{T}}

##
# We can identify a product space with equal types with a vector space of the same dimension and type
# The agreement is between NTuple{N,T} and SVector{N,T}.
##

embedding(::Type{ProductSpace{NTuple{N,S}}}, ::VectorSpace{N,T}) where {N,S,T} = promotes_to(S,T)
embedding(::VectorSpace{N,S}, ::Type{ProductSpace{NTuple{N,T}}}) where {N,S,T} = promotes_to(S,T)

function convert(::Type{ProductSpace{NTuple{N,T}}}, x::SVector{N,S}) where {N,T,S}
    NTuple{N,T}(x)
end

function convert(::Type{VectorSpace{N,T}}, x::NTuple{N,S}) where {N,T,S}
    SVector{N,T}(x)
end

# We prefer the SVector representation in case of mixing, since it acts like a vector.
isomorphism_promotion_rule(::Type{ProductSpace{NTuple{N,S}}}, ::VectorSpace{N,T}) where {N,S,T} = VectorSpace{N,T}
