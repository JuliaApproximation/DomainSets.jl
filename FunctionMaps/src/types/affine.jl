
"""
An affine map has the general form `y = A*x + b`.

We use `affinematrix(m)` and `affinevector(m)` to denote `A` and `b` respectively. Concrete
subtypes include linear maps of the form `y = A*x` and translations of the
form `y = x + b`.
"""
abstract type AbstractAffineMap{T} <: Map{T} end

unsafe_matrix(m::AbstractAffineMap) = m.A
unsafe_vector(m::AbstractAffineMap) = m.b

"Return the matrix `A` in the affine map `Ax+b`."
affinematrix(m::AbstractAffineMap) = to_matrix(domaintype(m), unsafe_matrix(m), unsafe_vector(m))

"Return the vector `b` in the affine map `Ax+b`."
affinevector(m::AbstractAffineMap) = to_vector(domaintype(m), unsafe_matrix(m), unsafe_vector(m))

applymap(m::AbstractAffineMap, x) = _affine_applymap(m, x, unsafe_matrix(m), unsafe_vector(m))
_affine_applymap(m, x, A, b) = A*x + b

applymap!(y, m::AbstractAffineMap, x) = _affine_applymap!(y, m, x, unsafe_matrix(m), unsafe_vector(m))
function _affine_applymap!(y, m, x, A, b)
    mul!(y, A, x)
    y .+= b
    y
end

isreal(m::AbstractAffineMap) = _affine_isreal(m, unsafe_matrix(m), unsafe_vector(m))
_affine_isreal(m, A, b) = isreal(A) && isreal(b)

jacobian(m::AbstractAffineMap{T}) where {T} = ConstantMap{T}(affinematrix(m))
jacobian(m::AbstractAffineMap, x) = affinematrix(m)

jacdet(m::AbstractAffineMap, x) = _affine_jacdet(m, x, unsafe_matrix(m))
_affine_jacdet(m, x, A) = det(A)
_affine_jacdet(m, x::Number, A::UniformScaling) = A.λ
_affine_jacdet(m, x::AbstractVector, A::Number) = A^length(x)
_affine_jacdet(m, x::AbstractVector, A::UniformScaling) = A.λ^length(x)

function diffvolume(m::AbstractAffineMap{T}) where T
    J = jacobian(m)
    c = sqrt(det(affinevector(J)'*affinevector(J)))
    ConstantMap{T}(c)
end

islinearmap(m::AbstractMap) = false
islinearmap(m::AbstractAffineMap) = _affine_islinearmap(m, unsafe_vector(m))
_affine_islinearmap(m, b) = all(b .== 0)

isaffinemap(m) = false
isaffinemap(m::Map) = islinearmap(m) || isconstantmap(m)
isaffinemap(m::AbstractAffineMap) = true

isequalmap(m1::AbstractAffineMap, m2::AbstractAffineMap) =
    affinematrix(m1) == affinematrix(m2) && affinevector(m1) == affinevector(m2)

isequalmap(m1::AbstractAffineMap, m2::IdentityMap) =
    islinearmap(m1) && affinematrix(m1) == affinematrix(m2)
isequalmap(m1::IdentityMap, m2::AbstractAffineMap) = isequalmap(m2, m1)

map_hash(m::AbstractAffineMap, h::UInt) = hashrec("AbstractAffineMap", affinematrix(m), affinevector(m), h)

mapsize(m::AbstractAffineMap) = _affine_mapsize(m, domaintype(m), unsafe_matrix(m), unsafe_vector(m))
_affine_mapsize(m, T, A::AbstractArray, b) = size(A)
_affine_mapsize(m, T, A::AbstractVector, b::AbstractVector) = (length(A),)
_affine_mapsize(m, T, A::Number, b::Number) = ()
_affine_mapsize(m, T, A::Number, b::AbstractVector) = (length(b),length(b))
_affine_mapsize(m, T, A::UniformScaling, b::Number) = ()
_affine_mapsize(m, T, A::UniformScaling, b) = (length(b),length(b))


Display.displaystencil(m::AbstractAffineMap) = vcat(["x -> "], map_stencil(m, 'x'))
show(io::IO, mime::MIME"text/plain", m::AbstractAffineMap) = composite_show(io, mime, m)

map_stencil(m::AbstractAffineMap, x) = _affine_map_stencil(m, x, unsafe_matrix(m), unsafe_vector(m))
_affine_map_stencil(m, x, A, b) = [A, " * ", x, " + ", b]
_affine_map_stencil(m, x, A, b::Real) =
    b >= 0 ? [A, " * ", x, " + ", b] :  [A, " * ", x, " - ", abs(b)]

map_stencil_broadcast(m::AbstractAffineMap, x) = _affine_map_stencil_broadcast(m, x, unsafe_matrix(m), unsafe_vector(m))
_affine_map_stencil_broadcast(m, x, A, b) = [A, " .* ", x, " .+ ", b]
_affine_map_stencil_broadcast(m, x, A::Number, b) = [A, " * ", x, " .+ ", b]

map_object_parentheses(m::AbstractAffineMap) = true
map_stencil_parentheses(m::AbstractAffineMap) = true

########################
# Linear maps: y = A*x
########################

"""
The supertype of all linear maps `y = A*x`.
Concrete subtypes may differ in how `A` is represented.
"""
abstract type LinearMap{T} <: AbstractAffineMap{T} end

mapsize(m::LinearMap) = _linearmap_size(m, domaintype(m), unsafe_matrix(m))
_linearmap_size(m, T, A) = size(A)
_linearmap_size(m, ::Type{T}, A::Number) where {T<:Number} = ()
_linearmap_size(m, ::Type{T}, A::AbstractVector) where {T<:Number} = (length(A),)
_linearmap_size(m, ::Type{T}, A::Number) where {N,T<:StaticVector{N}} = (N,N)

affinematrix(m::LinearMap) = to_matrix(domaintype(m), unsafe_matrix(m))
affinevector(m::LinearMap) = to_vector(domaintype(m), unsafe_matrix(m))

LinearMap(A::Number) = ScalarLinearMap(A)
LinearMap(A::SMatrix) = StaticLinearMap(A)
LinearMap(A::MMatrix) = StaticLinearMap(A)
LinearMap(A::Matrix) = VectorLinearMap(A)
LinearMap(A) = GenericLinearMap(A)

LinearMap{T}(A::Number) where {T<:Number} = _LinearMap(A, T, promote_type(T,typeof(A)))
_LinearMap(A::Number, ::Type{T}, ::Type{T}) where {T<:Number} = ScalarLinearMap{T}(A)
_LinearMap(A::Number, ::Type{T}, ::Type{S}) where {T<:Number,S} = GenericLinearMap{T}(A)
LinearMap{T}(A::SMatrix{M,N}) where {M,N,S,T <: SVector{N,S}} = StaticLinearMap{S}(A)
LinearMap{T}(A::MMatrix{M,N}) where {M,N,S,T <: SVector{N,S}} = StaticLinearMap{S}(A)
LinearMap{T}(A::Matrix) where {S,T <: Vector{S}} = VectorLinearMap{S}(A)
LinearMap{T}(A) where {T} = GenericLinearMap{T}(A)

# convenience functions
LinearMap(a::Number...) = LinearMap(promote(a...))
LinearMap(a::NTuple{N,T}) where {N,T <: Number} = LinearMap{SVector{N,T}}(Diagonal(SVector{N,T}(a)))

applymap(m::LinearMap, x) = _linear_applymap(m, x, unsafe_matrix(m))
_linear_applymap(m, x, A) = A*x

applymap!(y, m::LinearMap, x) = _linear_applymap!(y, m, x, unsafe_matrix(m))
_linear_applymap!(y, m, x, A) = mul!(y, A, x)

islinearmap(m::LinearMap) = true

isreal(m::LinearMap) = _linear_isreal(m, unsafe_matrix(m))
_linear_isreal(m, A) = isreal(A)

isequalmap(m1::LinearMap, m2::LinearMap) = affinematrix(m1) == affinematrix(m2)

# inverse should be called only on square maps, otherwise use
# leftinverse or rightinverse in order to use pinv instead of inv
inverse(m::LinearMap) = (@assert issquaremap(m); LinearMap(inv(m.A)))
inverse(m::LinearMap, x) = (@assert issquaremap(m); m.A \ x)

function leftinverse(m::LinearMap)
    @assert isoverdetermined(m)
    LinearMap(matrix_pinv(m.A))
end
function rightinverse(m::LinearMap)
    @assert isunderdetermined(m)
    LinearMap(matrix_pinv(m.A))
end
function leftinverse(m::LinearMap, x)
    @assert isoverdetermined(m)
    m.A \ x
end
function rightinverse(m::LinearMap, x)
    @assert isunderdetermined(m)
    m.A \ x
end

similarmap(m::LinearMap, ::Type{T}) where {T} = LinearMap{T}(m.A)

convert(::Type{Map}, a::Number) = LinearMap(a)
convert(::Type{Map{T}}, a::Number) where {T} = LinearMap{T}(a)

convert(::Type{LinearMap}, ::StaticIdentityMap{T}) where {T} = LinearMap{T}(1)
convert(::Type{LinearMap{T}}, ::StaticIdentityMap) where {T} = LinearMap{T}(1)
convert(::Type{AbstractAffineMap}, m::StaticIdentityMap) = convert(LinearMap, m)
convert(::Type{AbstractAffineMap{T}}, m::StaticIdentityMap) where {T} = convert(LinearMap, m)


map_stencil(m::LinearMap, x) = [unsafe_matrix(m), " * ", x]
map_stencil_broadcast(m::LinearMap, x) = _linear_map_stencil_broadcast(m, x, unsafe_matrix(m))
_linear_map_stencil_broadcast(m, x, A) = [A, " .* ", x]
_linear_map_stencil_broadcast(m, x, A::Number) = [A, " * ", x]


"A `GenericLinearMap` is a linear map `y = A*x` for any type of `A`."
struct GenericLinearMap{T,AA} <: LinearMap{T}
    A   ::  AA
end

"""
What is the suggested domaintype for a generic linear map `A*x` with
the given argument 'A'?
"""
glm_domaintype(A) = Any
glm_domaintype(A::Number) = typeof(A)
glm_domaintype(A::AbstractMatrix{T}) where T = Vector{T}
glm_domaintype(A::StaticMatrix{M,N,T}) where {M,N,T} = SVector{N,T}
glm_domaintype(A::AbstractVector{T}) where {T} = T
glm_domaintype(A::Diagonal{T,<:StaticVector{N,T}}) where {N,T} = SVector{N,T}

GenericLinearMap(A) = GenericLinearMap{glm_domaintype(A)}(A)

# Allow any A
GenericLinearMap{T}(A) where {T} = GenericLinearMap{T,typeof(A)}(A)

# Promote some eltypes if applicable
GenericLinearMap{T}(A::Number) where {T <: Number} =
    _GenericLinearMap(A, T, promote_type(T,typeof(A)))
_GenericLinearMap(A::Number, ::Type{T}, ::Type{T}) where {T <: Number} =
    GenericLinearMap{T,T}(A)
_GenericLinearMap(A::Number, ::Type{T}, ::Type{S}) where {S,T<:Number} =
    GenericLinearMap{T,typeof(A)}(A)

GenericLinearMap{T}(A::Number) where {T <: AbstractVector} =
    _GenericLinearMap(A, T, promote_type(eltype(T),typeof(A)))
_GenericLinearMap(A::Number, ::Type{T}, ::Type{S}) where {S,T<:AbstractVector{S}} =
    GenericLinearMap{T,S}(A)
_GenericLinearMap(A::Number, ::Type{T}, ::Type{S}) where {S,T<:AbstractVector} =
    GenericLinearMap{T,typeof(A)}(A)

GenericLinearMap{T}(A::AbstractMatrix{S}) where {S,T <: AbstractVector{S}} =
    GenericLinearMap{T,typeof(A)}(A)
GenericLinearMap{T}(A::AbstractMatrix{U}) where {S,T <: AbstractVector{S},U} =
    GenericLinearMap{T}(convert(AbstractMatrix{S}, A))
GenericLinearMap{T}(A::AbstractVector{T}) where {T} =
    GenericLinearMap{T,typeof(A)}(A)
GenericLinearMap{T}(A::AbstractVector{S}) where {S,T} =
    GenericLinearMap{T}(convert(AbstractVector{T}, A))
GenericLinearMap{T}(A::AbstractMatrix) where {T<:Number} =
    throw(ArgumentError("Linear map with matrix A can not have scalar domaintype."))

# Preserve the action on vectors with a number type
inverse(m::GenericLinearMap{T,AA}) where {T<:AbstractVector,AA<:Number} =
    LinearMap{T}(inv(m.A))
leftinverse(m::GenericLinearMap{T,AA}) where {T<:AbstractVector,AA<:Number} =
    LinearMap{T}(inv(m.A))
rightinverse(m::GenericLinearMap{T,AA}) where {T<:AbstractVector,AA<:Number} =
    LinearMap{T}(inv(m.A))

convert(::Type{Map}, A::UniformScaling) = GenericLinearMap{Vector{Any}}(A)
convert(::Type{Map{T}}, A::UniformScaling) where {T} = GenericLinearMap{T}(A)

"A `ScalarLinearMap` is a linear map `y = A*x` for scalars."
struct ScalarLinearMap{T} <: LinearMap{T}
    A   ::  T
end

ScalarLinearMap(A::Number) = ScalarLinearMap{typeof(A)}(A)

isreal(m::ScalarLinearMap{T}) where {T} = isrealtype(T)

show(io::IO, m::ScalarLinearMap) = show_scalar_linear_map(io, m.A)
show_scalar_linear_map(io, A::Real) = print(io, "x -> $(A) * x")
show_scalar_linear_map(io, A::Complex) = print(io, "x -> ($(A)) * x")


"A `VectorLinearMap` is a linear map `y = A*x` using vectors and matrices."
struct VectorLinearMap{T} <: LinearMap{Vector{T}}
    A   ::  Matrix{T}
end

VectorLinearMap(A::AbstractMatrix{T}) where {T} =
    VectorLinearMap{T}(A)


"A `StaticLinearMap` is a linear map `y = A*x` using static arrays."
struct StaticLinearMap{T,N,M,L} <: LinearMap{SVector{N,T}}
    A   ::  SMatrix{M,N,T,L}
end

StaticLinearMap(A::AbstractMatrix{T}) where {T} =
    StaticLinearMap{T}(A)

StaticLinearMap{T}(A::StaticMatrix{M,N,S}) where {M,N,T,S} =
    StaticLinearMap{T}(convert(AbstractMatrix{T}, A))
StaticLinearMap{T}(A::StaticMatrix{M,N,T}) where {M,N,T} =
    StaticLinearMap{T,M,N}(A)
StaticLinearMap{T,N,M}(A::AbstractMatrix) where {T,N,M} =
    StaticLinearMap{T,N,M,M*N}(A)

convert(::Type{Map{SVector{N,T}}}, m::VectorLinearMap) where {N,T} = StaticLinearMap{T,N,N}(m.A)


# Implement the interface for abstract arrays,
# representing the linear map x->A*x
MapStyle(A::AbstractArray) = IsMap()

Map(A::AbstractArray) = LinearMap(A)
Map{T}(A::AbstractArray) where T = LinearMap{T}(A)

domaintype(A::AbstractArray) = domaintype(Map(A))

applymap(A::AbstractArray, x) = A*x
mapsize(A::AbstractArray) = size(A)

islinearmap(A::AbstractArray) = true
isaffinemap(A::AbstractArray) = true
affinematrix(A::AbstractArray) = A
affinevector(A::AbstractArray) = zerovector(A)

inverse(A::AbstractMatrix) = inv(A)
inverse(A::AbstractMatrix, x) = A \ x

jacobian(A::AbstractMatrix) = ConstantMap{glm_domaintype(A)}(A)
jacobian(A::AbstractMatrix, x) = A

canonicalmap(A::AbstractArray) = Map(A)
canonicalmap(::Equal, A::AbstractArray) = Map(A)




##########################
# Translations: y = x + b
##########################


"A `Translation` represents the map `y = x + b`."
abstract type Translation{T} <: AbstractAffineMap{T} end

# unsafe_matrix(m::Translation) = I
affinematrix(m::Translation) = to_matrix(domaintype(m), LinearAlgebra.I, unsafe_vector(m))
affinevector(m::Translation) = unsafe_vector(m)

mapsize(m::Translation) = _translation_mapsize(m, domaintype(m), unsafe_vector(m))
_translation_mapsize(m, ::Type{T}, b::Number) where {T} = ()
_translation_mapsize(m, ::Type{T}, b) where {T} = (length(b),length(b))

map_stencil(m::Translation, x) = _translation_map_stencil(m, x, unsafe_vector(m))
_translation_map_stencil(m, x, b) = [x, " + ", b]
_translation_map_stencil(m, x, b::Real) =
    b >= 0 ? [x, " + ", b] : [x, " - ", abs(b)]
map_stencil_broadcast(m::Translation, x) = _translation_map_stencil_broadcast(m, x, unsafe_vector(m))
_translation_map_stencil_broadcast(m, x, b) = [x, " .+ ", b]

"Translation by a scalar value."
struct ScalarTranslation{T} <: Translation{T}
    b   ::  T
end

isreal(m::ScalarTranslation{T}) where {T} = isrealtype(T)

show(io::IO, m::ScalarTranslation) = show_scalar_translation(io, m.b)
show_scalar_translation(io, b::Real) = print(io, "x -> x", b < 0 ? " - " : " + ", abs(b))
show_scalar_translation(io, b) = print(io, "x -> x + ", b)


"Translation by a static vector."
struct StaticTranslation{T,N} <: Translation{SVector{N,T}}
    b   ::  SVector{N,T}
end

"Translation by a vector."
struct VectorTranslation{T} <: Translation{Vector{T}}
    b   ::  Vector{T}
end

"Translation by a generic vectorlike object."
struct GenericTranslation{T,B} <: Translation{T}
    b   ::  B
end

Translation(b::Number) = ScalarTranslation(b)
Translation(b::StaticVector) = StaticTranslation(b)
Translation(b::Vector) = VectorTranslation(b)
Translation(b) = GenericTranslation(b)

Translation{T}(b::Number) where {T<:Number} = ScalarTranslation{T}(b)
Translation{T}(b::AbstractVector) where {N,S,T<:StaticVector{N,S}} = StaticTranslation{S,N}(b)
Translation{T}(b::Vector) where {S,T<:Vector{S}} = VectorTranslation{S}(b)
Translation{T}(b) where {T} = GenericTranslation{T}(b)

jacdet(m::Translation{T}) where {T} = UnityMap{T,eltype(T)}()
jacdet(m::Translation{T}, x) where {T} = one(eltype(T))

isreal(m::Translation) = isreal(unsafe_vector(m))

isequalmap(m1::Translation, m2::Translation) = unsafe_vector(m1)==unsafe_vector(m2)

similarmap(m::Translation, ::Type{T}) where {T} = Translation{T}(m.b)

applymap(m::Translation, x) = _translation_applymap(m, x, unsafe_vector(m))
_translation_applymap(m, x, b) = x + b
applymap!(y, m::Translation, x) = _translation_applymap!(y, m, x, unsafe_vector(m))
_translation_applymap!(y, m, x, b) = y .= x .+ m.b

inverse(m::Translation{T}) where {T} = Translation{T}(-m.b)
inverse(m::Translation, x) = x - m.b

ScalarTranslation(b::Number) = ScalarTranslation{typeof(b)}(b)

StaticTranslation(b::AbstractVector{T}) where {T} = StaticTranslation{T}(b)

StaticTranslation{T}(b::StaticVector{N}) where {N,T} =
    StaticTranslation{T,N}(b)

VectorTranslation(b::AbstractVector{T}) where {T} = VectorTranslation{T}(b)

GenericTranslation(b) = GenericTranslation{typeof(b)}(b)
GenericTranslation(b::AbstractVector{T}) where {T} =
    GenericTranslation{Vector{T}}(b)
GenericTranslation(b::StaticVector{N,T}) where {N,T} =
    GenericTranslation{SVector{N,T}}(b)

GenericTranslation{T}(b) where {T} = GenericTranslation{T,typeof(b)}(b)
GenericTranslation{T}(b::Number) where {T<:Number} =
    GenericTranslation{T,T}(b)


############################
# Affine maps: y = A*x + b
############################


"""
The supertype of all affine maps that store `A` and `b`.
Concrete subtypes differ in how `A` and `b` are represented.
"""
abstract type AffineMap{T} <: AbstractAffineMap{T} end

AffineMap(A::Number, b::Number) = ScalarAffineMap(A, b)
AffineMap(A::StaticMatrix, b::StaticVector) = StaticAffineMap(A, b)
AffineMap(A::Matrix, b::Vector) = VectorAffineMap(A, b)
AffineMap(A::UniformScaling{Bool}, b::Number) = ScalarAffineMap(one(b), b)
AffineMap(A, b) = GenericAffineMap(A, b)

AffineMap{T}(A::Number, b::Number) where {T<:Number} = ScalarAffineMap{T}(A, b)
AffineMap{T}(A::AbstractMatrix, b::AbstractVector) where {N,S,T<:SVector{N,S}} = StaticAffineMap{S,N}(A, b)
AffineMap{T}(A::Matrix, b::Vector) where {S,T<:Vector{S}} = VectorAffineMap{S}(A, b)
AffineMap{T}(A::UniformScaling{Bool}, b::Number) where {T} = ScalarAffineMap{T}(one(T), b)
AffineMap{T}(A, b) where {T} = GenericAffineMap{T}(A, b)

similarmap(m::AffineMap, ::Type{T}) where {T} = AffineMap{T}(m.A, m.b)

convert(::Type{AffineMap}, m) = (@assert isaffinemap(m); AffineMap(affinematrix(m), affinevector(m)))
convert(::Type{AffineMap{T}}, m) where {T} = (@assert isaffinemap(m); AffineMap{T}(affinematrix(m), affinevector(m)))
# avoid ambiguity errors with convert(::Type{T}, x::T) in Base:
convert(::Type{AffineMap}, m::AffineMap) = m
convert(::Type{AffineMap{T}}, m::AffineMap{T}) where T = m

# If y = A*x+b, then x = inv(A)*(y-b) = inv(A)*y - inv(A)*b
inverse(m::AffineMap) = (@assert issquaremap(m); AffineMap(inv(m.A), -inv(m.A)*m.b))
inverse(m::AffineMap, x) = (@assert issquaremap(m); m.A \ (x-m.b))

function leftinverse(m::AffineMap)
    @assert isoverdetermined(m)
    pA = matrix_pinv(m.A)
    AffineMap(pA, -pA*m.b)
end
function rightinverse(m::AffineMap)
    @assert isunderdetermined(m)
    pA = matrix_pinv(m.A)
    AffineMap(pA, -pA*m.b)
end
function leftinverse(m::AffineMap, x)
    @assert isoverdetermined(m)
    m.A \ (x-m.b)
end
function rightinverse(m::AffineMap, x)
    @assert isunderdetermined(m)
    m.A \ (x-m.b)
end


"An affine map for any combination of types of `A` and `b`."
struct GenericAffineMap{T,AA,B} <: AffineMap{T}
    A   ::  AA
    b   ::  B
end

GenericAffineMap(A, b) = GenericAffineMap{typeof(b)}(A, b)
GenericAffineMap(A::AbstractVector{S}, b::AbstractVector{T}) where {S,T} =
    GenericAffineMap{promote_type(S,T)}(A, b)
GenericAffineMap(A::AbstractArray{S}, b::AbstractVector{T}) where {S,T} =
    GenericAffineMap{Vector{promote_type(S,T)}}(A, b)
GenericAffineMap(A::StaticMatrix{M,N,S}, b::StaticVector{M,T}) where {M,N,S,T} =
    GenericAffineMap{SVector{N,promote_type(S,T)}}(A, b)
GenericAffineMap(A::StaticMatrix{M,N,S}, b::AbstractVector{T}) where {M,N,S,T} =
    GenericAffineMap{SVector{N,promote_type(S,T)}}(A, b)
GenericAffineMap(A::S, b::AbstractVector{T}) where {S<:Number,T} =
    GenericAffineMap{Vector{promote_type(S,T)}}(A, b)
GenericAffineMap(A::S, b::StaticVector{N,T}) where {S<:Number,N,T} =
    GenericAffineMap{SVector{N,promote_type(S,T)}}(A, b)
GenericAffineMap(A::UniformScaling{Bool}, b) =
    GenericAffineMap(UniformScaling{eltype(b)}(1), b)


# Fallback routine for generic A and b, special cases follow
GenericAffineMap{T}(A, b) where {T} = GenericAffineMap{T,typeof(A),typeof(b)}(A, b)

GenericAffineMap{T}(A::AbstractVector{S}, b::AbstractVector{U}) where {T<:Number,S,U} =
    GenericAffineMap{T}(convert(AbstractVector{T}, A), convert(AbstractVector{T}, b))
GenericAffineMap{T}(A::AbstractVector{T}, b::AbstractVector{T}) where {T<:Number} =
    GenericAffineMap{T,typeof(A),typeof(b)}(A, b)
GenericAffineMap{T}(A::Number, b) where {T} = GenericAffineMap{T,eltype(T),typeof(b)}(A, b)
GenericAffineMap{T}(A::Number, b::AbstractVector) where {N,S,T <: StaticVector{N,S}} =
    GenericAffineMap{T,S,SVector{N,S}}(A, b)
# Promote element types of abstract arrays
GenericAffineMap{T}(A::AbstractArray, b::AbstractVector) where {S,T<:AbstractVector{S}} =
    GenericAffineMap{T}(convert(AbstractArray{eltype(T)},A), convert(AbstractVector{eltype(T)}, b))
GenericAffineMap{T}(A::AbstractArray{S}, b::AbstractVector{S}) where {S,T<:AbstractVector{S}} =
    GenericAffineMap{T,typeof(A),typeof(b)}(A, b)
GenericAffineMap{T}(A::UniformScaling{Bool}, b::AbstractVector) where {S,T<:AbstractVector{S}} =
    GenericAffineMap{T}(A*one(S), convert(AbstractVector{S}, b))
GenericAffineMap{T}(A::UniformScaling{S}, b::AbstractVector{S}) where {S,T<:AbstractVector{S}} =
    GenericAffineMap{T,typeof(A),typeof(b)}(A, b)


similarmap(m::GenericAffineMap, ::Type{T}) where {T} = AffineMap{T}(m.A, m.b)

convert(::Type{GenericAffineMap{T}}, m::GenericAffineMap) where {T} =
    GenericAffineMap{T}(m.A, m.b)



"An affine map with scalar representation."
struct ScalarAffineMap{T} <: AffineMap{T}
    A   ::  T
    b   ::  T
end

ScalarAffineMap(A, b) = ScalarAffineMap(promote(A, b)...)

isreal(m::ScalarAffineMap{T}) where {T} = isrealtype(T)

show(io::IO, m::ScalarAffineMap) = show_scalar_affine_map(io, m.A, m.b)
show_scalar_affine_map(io, A::Real, b::Real) = print(io, "x -> $(A) * x", b < 0 ? " - " : " + ", abs(b))
show_scalar_affine_map(io, A::Complex, b::Complex) = print(io, "x -> ($(A)) * x + ", b)
show_scalar_affine_map(io, A, b) = print(io, "x -> ($(A)) * x + $(b)")


convert(::Type{ScalarAffineMap{T}}, m::ScalarAffineMap) where {T} =
    ScalarAffineMap{T}(m.A, m.b)

"An affine map with array and vector representation."
struct VectorAffineMap{T} <: AffineMap{Vector{T}}
    A   ::  Matrix{T}
    b   ::  Vector{T}
end

VectorAffineMap(A::AbstractArray{T}, b::AbstractVector{T}) where {T} =
    VectorAffineMap{T}(A, b)
function VectorAffineMap(A::AbstractArray{S}, b::AbstractVector{T}) where {S,T}
    U = promote_type(S,T)
    VectorAffineMap(convert(AbstractArray{U}, A), convert(AbstractVector{U}, b))
end

convert(::Type{VectorAffineMap{T}}, m::VectorAffineMap) where {T} =
    VectorAffineMap{T}(m.A, m.b)



"An affine map with representation using static arrays."
struct StaticAffineMap{T,N,M,L} <: AffineMap{SVector{N,T}}
    A   ::  SMatrix{M,N,T,L}
    b   ::  SVector{M,T}
end

# Constructors:
# - first, we deduce T
StaticAffineMap(A::AbstractMatrix{T}, b::AbstractVector{T}) where {T} =
    StaticAffineMap{T}(A, b)
function StaticAffineMap(A::AbstractMatrix{S}, b::AbstractVector{T}) where {S,T}
    U = promote_type(S,T)
    StaticAffineMap(convert(AbstractMatrix{U}, A), convert(AbstractVector{U}, b))
end

StaticAffineMap{T}(A::AbstractMatrix, b::AbstractVector) where {T} =
    StaticAffineMap{T}(convert(AbstractMatrix{T}, A), convert(AbstractVector{T}, b))

# - then, we determine N and/or M, from the arguments
function StaticAffineMap{T}(A::AbstractMatrix{T}, b::StaticVector{M,T}) where {T,M}
    @assert size(A) == (M,M)
    StaticAffineMap{T,M,M}(A, b)
end
StaticAffineMap{T}(A::StaticMatrix{M,N,T}, b::AbstractVector) where {T,N,M} =
    StaticAffineMap{T,N,M}(A, b)
StaticAffineMap{T}(A::StaticMatrix{M,N,T}, b::StaticVector{M,T}) where {T,N,M} =
    StaticAffineMap{T,N,M}(A, b)
# line below catches ambiguity error
StaticAffineMap{T}(A::StaticMatrix{M1,N,T}, b::StaticVector{M2,T}) where {T,N,M1,M2} =
    throw(ArgumentError("Non-matching dimensions"))
StaticAffineMap{T,N}(A::AbstractMatrix, b::AbstractVector) where {T,N} =
    StaticAffineMap{T,N,N}(A, b)
StaticAffineMap{T,N}(A::StaticMatrix{M,N}, b::AbstractVector) where {T,N,M} =
    StaticAffineMap{T,N,M}(A, b)

# - finally invoke the constructor (and implicitly convert the data if necessary)
StaticAffineMap{T,N,M}(A::AbstractMatrix, b::AbstractVector) where {T,N,M} =
    StaticAffineMap{T,N,M,M*N}(A, b)

convert(::Type{Map{SVector{N,T}}}, m::VectorAffineMap) where {N,T} =
    StaticAffineMap{T,N}(m.A, m.b)
