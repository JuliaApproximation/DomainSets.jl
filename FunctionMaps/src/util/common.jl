
isrealtype(::Type{<:Real}) = true
isrealtype(::Type{<:Complex}) = false
isrealtype(::Type{T}) where {T} = isrealtype(eltype(T))

const StaticTypes = Union{Number,<:StaticVector{N} where N,<:NTuple{N,Any} where N}

"Apply the `hash` function recursively to the given arguments."
hashrec(x) = hash(x)
hashrec(x, args...) = hash(x, hashrec(args...))

# Workaround for #88, manually compute the hash of an array using all its elements
hashrec() = zero(UInt)
function hashrec(A::AbstractArray, args...)
	h = hash(size(A))
	for x in A
		h = hash(x, h)
	end
	hash(h, hashrec(args...))
end

"What is the euclidean dimension of the given type (if applicable)?"
euclideandimension(::Type{T}) where {T <: Number} = 1
euclideandimension(::Type{T}) where {N,T <: StaticVector{N}} = N
euclideandimension(::Type{T}) where {N,T <: NTuple{N,Any}} = N
euclideandimension(::Type{T}) where {T} =
	throw(ArgumentError("Don't know the euclidean dimension of $(T)."))

#################
# Element type
#################

"""
	convert_eltype(T, x)

Convert `x` such that its `eltype` equals `T`.
"""
convert_eltype(::Type{T}, d::AbstractArray) where {T} = convert(AbstractArray{T}, d)
convert_eltype(::Type{T}, d::AbstractRange) where {T} = map(T, d)
convert_eltype(::Type{T}, d::Set) where {T} = convert(Set{T}, d)
convert_eltype(::Type{T}, d::Number) where {T} = convert(T, d)

"Are the given element types promotable to a concrete supertype?"
promotable_eltypes(types...) = isconcretetype(promote_type(types...))
promotable_eltypes(::Type{S}, ::Type{T}) where {S<:AbstractVector,T<:AbstractVector} =
    promotable_eltypes(eltype(S), eltype(T))

#################
# Precision type
#################

"""
	prectype(x[, ...])

The floating point precision type associated with the argument(s).
"""
prectype(x) = prectype(typeof(x))
prectype(::Type{T}) where T = Any		# fallback
prectype(::Type{T}) where {T<:AbstractFloat} = T
prectype(::Type{T}) where {T<:Number} = prectype(float(T))
prectype(::Type{<:AbstractArray{T}}) where T = prectype(T)
prectype(::Type{<:Complex{T}}) where {T} = prectype(T)
prectype(::Type{NTuple{N,T}}) where {N,T} = prectype(T)
prectype(::Type{Tuple{A}}) where {A} = prectype(A)
prectype(::Type{Tuple{A,B}}) where {A,B} = prectype(A,B)
prectype(::Type{Tuple{A,B,C}}) where {A,B,C} = prectype(A,B,C)
@generated function prectype(T::Type{<:Tuple})
    quote $(promote_type(map(prectype, T.parameters[1].parameters)...)) end
end

prectype(a...) = promote_type(map(prectype, a)...)

"""
	convert_prectype(T, x)

Convert `x` such that its `prectype` equals `T`.
"""
convert_prectype(::Type{T}, x) where {T} =
	convert_eltype(to_prectype(T, eltype(x)), x)

"""
	to_prectype(U, T)

Return the type to which `T` can be converted, such that the `prectype` becomes `U`.
"""
to_prectype(::Type{U}, ::Type{T}) where {T,U} = throw(ArgumentError("Don't know how to convert the prectype of $(T) to $(U)."))
to_prectype(::Type{U}, ::Type{T}) where {T <: Real,U <: Real} = U
to_prectype(::Type{U}, ::Type{Complex{T}}) where {T <: Real,U <: Real} = Complex{U}
to_prectype(::Type{U}, ::Type{Vector{T}}) where {T,U} = Vector{to_prectype(U,T)}
to_prectype(::Type{U}, ::Type{<:StaticVector{N,T}}) where {N,T,U} = SVector{N,to_prectype(U,T)}
to_prectype(::Type{U}, ::Type{MVector{N,T}}) where {N,T,U} = MVector{N,to_prectype(U,T)}
to_prectype(::Type{U}, ::Type{<:StaticMatrix{M,N,T}}) where {M,N,T,U} = SMatrix{M,N,to_prectype(U,T)}
to_prectype(::Type{U}, ::Type{MMatrix{M,N,T}}) where {M,N,T,U} = MMatrix{M,N,to_prectype(U,T)}

"""
	promote_prectype(a, b[, ...])

Promote the precision types of the arguments to a joined supertype.
"""
promote_prectype(a) = a
promote_prectype(a, b) = _promote_prectype(prectype(a,b), a, b)
promote_prectype(a, b, c...) = _promote_prectype(prectype(a,b,c...), a, b, c...)
_promote_prectype(U, a) = convert_prectype(U, a)
_promote_prectype(U, a, b) = convert_prectype(U, a), convert_prectype(U, b)
_promote_prectype(U, a, b, c...) =
    (convert_prectype(U, a), convert_prectype(U, b), _promote_prectype(U, c...)...)

#################
# Numeric type
#################

"The numeric element type of `x` in a Euclidean space."
numtype(x) = numtype(typeof(x))
numtype(::Type{T}) where T = Any
numtype(::Type{T}) where {T<:Number} = T
numtype(::Type{<:AbstractArray{T}}) where T = T
numtype(::Type{NTuple{N,T}}) where {N,T} = T
numtype(::Type{Tuple{A,B}}) where {A,B} = promote_type(numtype(A), numtype(B))
numtype(::Type{Tuple{A,B,C}}) where {A,B,C} = promote_type(numtype(A), numtype(B), numtype(C))
numtype(::Type{Tuple{A,B,C,D}}) where {A,B,C,D} = promote_type(numtype(A), numtype(B), numtype(C), numtype(D))
@generated function numtype(T::Type{<:Tuple})
    quote $(promote_type(T.parameters[1].parameters...)) end
end

numtype(a...) = promote_type(map(numtype, a)...)

"""
	convert_numtype(T, x)

Convert `x` such that its `numtype` equals `T`.
"""
convert_numtype(::Type{T}, x) where {T} =
	convert_eltype(to_numtype(T, eltype(x)), x)

"""
	to_numtype(U, T)

Return the type to which `T` can be converted, such that the `numtype` becomes `U`.
"""
to_numtype(::Type{U}, ::Type{T}) where {T,U} = throw(ArgumentError("Don't know how to convert the numtype of $(T) to $(U)."))
to_numtype(::Type{U}, ::Type{T}) where {T <: Number,U <: Number} = U
to_numtype(::Type{U}, ::Type{Vector{T}}) where {T,U} = Vector{to_numtype(U,T)}
to_numtype(::Type{U}, ::Type{<:StaticVector{N,T}}) where {N,T,U} = SVector{N,to_numtype(U,T)}
to_numtype(::Type{U}, ::Type{MVector{N,T}}) where {N,T,U} = MVector{N,to_numtype(U,T)}
to_numtype(::Type{U}, ::Type{<:StaticMatrix{M,N,T}}) where {M,N,T,U} = SMatrix{M,N,to_numtype(U,T)}
to_numtype(::Type{U}, ::Type{MMatrix{M,N,T}}) where {M,N,T,U} = MMatrix{M,N,to_numtype(U,T)}

"""
	promote_numtype(a, b[, ...])

Promote the numeric types of the arguments to a joined supertype.
"""
promote_numtype(a) = a
promote_numtype(a, b) = _promote_numtype(numtype(a,b), a, b)
promote_numtype(a, b, c...) = _promote_numtype(numtype(a,b,c...), a, b, c...)
_promote_numtype(U, a) = convert_numtype(U, a)
_promote_numtype(U, a, b) = convert_numtype(U, a), convert_numtype(U, b)
_promote_numtype(U, a, b, c...) =
    (convert_numtype(U, a), convert_numtype(U, b), _promote_numtype(U, c...)...)

## Conversion from nested vectors to flat vectors and back

"""
Convert a vector from a cartesian format to a nested tuple according to the
given dimensions.

For example:
`convert_fromcartesian([1,2,3,4,5], Val{(2,2,1)}()) -> ([1,2],[3,4],5)`
"""
@generated function convert_fromcartesian(x::AbstractVector, ::Val{DIM}) where {DIM}
	dimsum = [0; cumsum([d for d in DIM])]
	E = Expr(:tuple, [ (dimsum[i+1]-dimsum[i] > 1 ? Expr(:call, :SVector, [:(x[$j]) for j = dimsum[i]+1:dimsum[i+1]]...) : :(x[$(dimsum[i+1])])) for i in 1:length(DIM)]...)
	return quote $(E) end
end

"The inverse function of `convert_fromcartesian`."
@generated function convert_tocartesian(x, ::Val{DIM}) where {DIM}
    dimsum = [0; cumsum([d for d in DIM])]
    E = vcat([[:(x[$i][$j]) for j in 1:DIM[i]] for i in 1:length(DIM)]...)
    quote SVector($(E...)) end
end

# we use matrix_pinv rather than pinv to preserve static matrices
matrix_pinv(A) = pinv(A)
# matrix_pinv(A::SMatrix{M,N}) where {M,N} = SMatrix{N,M}(pinv(A))
matrix_pinv(A::StaticVector{N}) where {N} = convert(Transpose{Float64, SVector{N,Float64}}, pinv(A))

## Composite objects

"Factors of a product-like composite object (equivalent to `components(d)`)."
function factors end
"The number of factors of a product-like composite object."
nfactors(d) = length(factors(d))
"Factor `I...` of a product-like composite object."
factor(d, I...) = getindex(factors(d), I...)
