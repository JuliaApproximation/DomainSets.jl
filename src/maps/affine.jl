
import Base: convert

# These are temporary definitions until StaticArrays updates
Base.convert(::Type{AbstractArray{T}}, v::SVector{N,S}) where {N,T,S} =
    convert(SVector{N,T}, v)
Base.convert(::Type{AbstractVector{T}}, v::SVector{N,S}) where {N,T,S} =
    convert(SVector{N,T}, v)
Base.convert(::Type{AbstractArray{T}}, v::SMatrix{M,N,S}) where {M,N,T,S} =
    convert(SMatrix{M,N,T}, v)
Base.convert(::Type{AbstractMatrix{T}}, v::SMatrix{M,N,S}) where {M,N,T,S} =
    convert(SMatrix{M,N,T}, v)

tomatrix(::Type{T}, a::Number) where {T<:Number} = a
tomatrix(::Type{SVector{N,T}}, a::Number) where {N,T} = a * one(SMatrix{N,N,T})
tomatrix(::Type{<:AbstractVector{T}}, a::Number) where {T} = a * I

"""
An affine map has the general form `y = A*x + b`.

We use `matrix(m)` and `vector(m)` to denote `A` and `b` respectively. Concrete
subtypes include linear maps of the form `y = A*x` and translations of the
form `y = x + b`.
"""
abstract type AbstractAffineMap{T} <: Map{T} end

"Return the matrix `A` in the affine map `Ax+b`."
matrix(m::AbstractAffineMap) = copy(unsafe_matrix(m))

"Return the vector `b` in the affine map `Ax+b`."
vector(m::AbstractAffineMap) = copy(unsafe_vector(m))

applymap(m::AbstractAffineMap, x) = _applymap(m, x, matrix(m), vector(m))
_applymap(m::AbstractAffineMap, x, A, b) = A*x+b

applymap!(y, m::AbstractAffineMap, x) = _applymap!(y, m, x, matrix(m), vector(m))
function applymap!(y, m::AbstractAffineMap, x, A, b)
    mul!(y, A, x)
    y .+= b
    y
end

isreal(m::AbstractAffineMap) = _isreal(m, matrix(m), vector(m))
_isreal(m::AbstractAffineMap, A) = isreal(A)
_isreal(m::AbstractAffineMap, A, b) = isreal(A) && isreal(b)

jacobian(m::AbstractAffineMap{T}) where {T} = ConstantMap{T}(matrix(m))
jacobian(m::AbstractAffineMap, x) = matrix(m)

jacdet(m::AbstractAffineMap, x) = _jacdet(m, x, matrix(m))
_jacdet(m::AbstractAffineMap, x, A) = det(A)
_jacdet(m::AbstractAffineMap, x::AbstractVector, A::Number) = det(A)^length(x)

islinear(m::AbstractMap) = false
islinear(m::AbstractAffineMap) = _islinear(m, vector(m))
_islinear(m::AbstractAffineMap, b) = all(b .== 0)

isaffine(m::AbstractMap) = islinear(m) || isconstant(m)
isaffine(m::AbstractAffineMap) = true

==(m1::AbstractAffineMap, m2::AbstractAffineMap) =
    (unsafe_matrix(m1) == unsafe_matrix(m2)) && (unsafe_vector(m1)==unsafe_vector(m2))

size(m::AbstractAffineMap) = _size(m, matrix(m), vector(m))
_size(m, A, b) = size(A)
_size(m, A::Number, b::Number) = (1,1)
_size(m, A::Number, b::AbstractVector) = (length(b),length(b))

unsafe_matrix(m::AbstractAffineMap) = m.A
unsafe_vector(m::AbstractAffineMap) = m.b


########################
# Linear maps: y = A*x
########################


"""
The supertype of all linear maps `y = A*x`.
Concrete subtypes may differ in how `A` is represented.
"""
abstract type LinearMap{T} <: AbstractAffineMap{T} end

vector(m::LinearMap{T}) where {T<:SVector} = zero(T)
vector(m::LinearMap{T}) where {T<:Number} = zero(T)
vector(m::LinearMap{T}) where {T<:AbstractVector} = zeros(eltype(T),size(m.A,1))

LinearMap(A::Number) = ScalarLinearMap(A)
LinearMap(A::SMatrix) = StaticLinearMap(A)
LinearMap(A::Matrix) = VectorLinearMap(A)
LinearMap(A) = GenericLinearMap(A)

LinearMap{T}(A::Number) where {T<:Number} = ScalarLinearMap{T}(A)
LinearMap{T}(A::SMatrix{M,N}) where {M,N,S,T <: SVector{N,S}} = StaticLinearMap{S}(A)
LinearMap{T}(A::Matrix) where {S,T <: Vector{S}} = VectorLinearMap{S}(A)
LinearMap{T}(A) where {T} = GenericLinearMap{T}(A)

applymap(m::LinearMap, x) = m.A * x
applymap!(y, m::LinearMap, x) = mul!(y, m.A, x)

islinear(m::LinearMap) = true

isreal(m::LinearMap) = _isreal(m, unsafe_matrix(m))

size(m::LinearMap) = _size(m, unsafe_matrix(m), 0)

==(m1::LinearMap, m2::LinearMap) = unsafe_matrix(m1) == unsafe_matrix(m2)

jacdet(m::LinearMap, x) = _jacdet(m, x, unsafe_matrix(m))

inv(m::LinearMap) = LinearMap(inv(m.A))
inverse(m::LinearMap, x) = m.A \ x

function leftinverse(m::LinearMap)
    M, N = size(m)
    M < N && error("No left inverse exists for $(m)")
    LinearMap(pinv(m.A))
end
function rightinverse(m::LinearMap)
    M, N = size(m)
    M > N && error("No right inverse exists for $(m)")
    LinearMap(pinv(m.A))
end
function leftinverse(m::LinearMap, x)
    M, N = size(m)
    @assert M >= N
    m.A \ x
end
function rightinverse(m::LinearMap, x)
    M, N = size(m)
    @assert M <= N
    m.A \ x
end

convert(::Type{Map{T}}, m::LinearMap{T}) where {T} = m
convert(::Type{Map{T}}, m::LinearMap{S}) where {S,T} = LinearMap{T}(m.A)

convert(::Type{Map{T}}, a::Number) where {T} = LinearMap{T}(a)

convert(::Type{AbstractAffineMap{T}}, ::IdentityMap) where {T} = LinearMap{T}(1)
convert(::Type{LinearMap{T}}, ::IdentityMap) where {T} = LinearMap{T}(1)


"A `GenericLinearMap` is a linear map `y = A*x` for any type of `A`."
struct GenericLinearMap{T,AA} <: LinearMap{T}
    A   ::  AA
end

GenericLinearMap(A::Number) = GenericLinearMap{typeof(A)}(A)
GenericLinearMap(A::AbstractMatrix{T}) where {T} =
    GenericLinearMap{Vector{T}}(A)
GenericLinearMap(A::SMatrix{M,N,T}) where {M,N,T} =
    GenericLinearMap{SVector{N,T}}(A)

# Allow any A
GenericLinearMap{T}(A) where {T} = GenericLinearMap{T,typeof(A)}(A)

GenericLinearMap{T}(A::Number) where {T <: Number} = GenericLinearMap{T,T}(A)
GenericLinearMap{T}(A::Number) where {S,T <: AbstractVector{S}} = GenericLinearMap{T,S}(A)
GenericLinearMap{T}(A::AbstractMatrix{S}) where {S,T <: AbstractVector{S}} =
    GenericLinearMap{T,typeof(A)}(A)
GenericLinearMap{T}(A::AbstractMatrix{U}) where {S,T <: AbstractVector{S},U} =
    GenericLinearMap{T}(convert(AbstractMatrix{S}, A))

# Preserve the action on vectors with a number type
inv(m::GenericLinearMap{T,AA}) where {T<:AbstractVector,AA<:Number} = LinearMap{T}(inv(m.A))


"A `ScalarLinearMap` is a linear map `y = A*x` for scalars."
struct ScalarLinearMap{T} <: LinearMap{T}
    A   ::  T
end

ScalarLinearMap(A::Number) = ScalarLinearMap{typeof(A)}(A)



"A `VectorLinearMap` is a linear map `y = A*x` using vectors and arrays."
struct VectorLinearMap{T} <: LinearMap{Vector{T}}
    A   ::  Array{T}
end

VectorLinearMap(A::AbstractMatrix{T}) where {T} =
    VectorLinearMap{T}(A)


"A `StaticLinearMap` is a linear map `y = A*x` using static arrays."
struct StaticLinearMap{T,N,M,L} <: LinearMap{SVector{N,T}}
    A   ::  SMatrix{M,N,T,L}
end

StaticLinearMap(A::AbstractMatrix{T}) where {T} =
    StaticLinearMap{T}(A)

StaticLinearMap{T}(A::SMatrix{M,N,S}) where {M,N,T,S} =
    StaticLinearMap{T}(convert(AbstractMatrix{T}, A))
StaticLinearMap{T}(A::SMatrix{M,N,T}) where {M,N,T} =
    StaticLinearMap{T,M,N}(A)
StaticLinearMap{T,N,M}(A::AbstractMatrix) where {T,N,M} =
    StaticLinearMap{T,N,M,M*N}(A)

convert(::Type{Map{SVector{N,T}}}, m::VectorLinearMap) where {N,T} = StaticLinearMap{T,N,N}(m.A)


##########################
# Translations: y = x + b
##########################


"A `Translation` represents the map `y = x + b`."
abstract type Translation{T} <: AbstractAffineMap{T} end

"Translation by a scalar value."
struct ScalarTranslation{T} <: Translation{T}
    b   ::  T
end

"Translation by a static vector."
struct StaticTranslation{T,N} <: Translation{SVector{N,T}}
    b   ::  SVector{N,T}
end

"Translation by a vector."
struct VectorTranslation{T} <: Translation{Vector{T}}
    b   ::  Vector{T}
end

"Translation by a generic vectorlike object."
struct GenericTranslation{T,B} <: Translation{Vector{T}}
    b   ::  B
end

Translation(b::Number) = ScalarTranslation(b)
Translation(b::SVector) = StaticTranslation(b)
Translation(b::Vector) = VectorTranslation(b)
Translation(b) = GenericTranslation(b)

Translation{T}(b::Number) where {T<:Number} = ScalarTranslation{T}(b)
Translation{T}(b::AbstractVector) where {N,S,T<:SVector{N,S}} = StaticTranslation{S,N}(b)
Translation{T}(b::Vector) where {S,T<:Vector{S}} = VectorTranslation{S}(b)
Translation{T}(b) where {T} = GenericTranslation{T}(b)

size(m::Translation) = (length(unsafe_vector(m)), length(unsafe_vector(m)))

matrix(m::Translation) = identitymatrix(m)

jacdet(m::Translation, x) = one(prectype(m))

isreal(m::Translation) = isreal(unsafe_vector(m))

==(m1::Translation, m2::Translation) = unsafe_vector(m1)==unsafe_vector(m2)


convert(::Type{Map{T}}, m::Translation{T}) where {T} = m
convert(::Type{Map{T}}, m::Translation{S}) where {S,T} = Translation{T}(m.b)

convert(::Type{Map{SVector{N,T}}}, m::Translation) where {N,T} = StaticTranslation{T,N}(m.b)

applymap(m::Translation, x) = x + m.b
applymap!(y, m, x) = y .+= x .+ m.b

inv(m::Translation{T}) where {T} = Translation{T}(-m.b)

inverse(m::Translation, x) = x - m.b


ScalarTranslation(b::Number) = ScalarTranslation{typeof(b)}(b)

StaticTranslation(b::AbstractVector{T}) where {T} = StaticTranslation{T}(b)

StaticTranslation{T}(b::AbstractVector{S}) where {S,T} =
    StaticTranslation{T}(convert(AbstractVector{T}, b))
StaticTranslation{T}(b::SVector{N,T}) where {N,T} =
    StaticTranslation{T,N}(b)

VectorTranslation(b::AbstractVector{T}) where {T} = VectorTranslation{T}(b)

GenericTranslation(b) = GenericTranslation{typeof(b)}(b)

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

applymap(m::AffineMap, x) = _applymap(m, x, unsafe_matrix(m), unsafe_vector(m))
applymap!(y, m::AffineMap, x) = _applymap!(y, m, x, unsafe_matrix(m), unsafe_vector(m))

AffineMap(A::Number, b::Number) = ScalarAffineMap(A, b)
AffineMap(A::SMatrix, b::SVector) = StaticAffineMap(A, b)
AffineMap(A::Matrix, b::Vector) = VectorAffineMap(A, b)
AffineMap(A, b) = GenericAffineMap(A, b)

AffineMap{T}(A::Number, b::Number) where {T<:Number} = ScalarAffineMap{T}(A, b)
AffineMap{T}(A::AbstractMatrix, b::AbstractVector) where {N,S,T<:SVector{N,S}} = StaticAffineMap{S,N}(A, b)
AffineMap{T}(A::Matrix, b::Vector) where {S,T<:Vector{S}} = VectorAffineMap{S}(A, b)
AffineMap{T}(A, b) where {T} = GenericAffineMap{T}(A, b)

convert(::Type{Map{T}}, m::AffineMap{T}) where {T} = m
convert(::Type{Map{T}}, m::AffineMap{S}) where {S,T} = AffineMap{T}(m.A, m.b)

isreal(m::AffineMap) = _isreal(m, unsafe_matrix(m), unsafe_vector(m))
islinear(m::AffineMap) = _islinear(m, unsafe_vector(m))

==(m1::AffineMap, m2::AffineMap) =
    (unsafe_matrix(m1) == unsafe_matrix(m2)) && (unsafe_vector(m1)==unsafe_vector(m2))

size(m::AffineMap) = _size(m, unsafe_matrix(m), unsafe_vector(m))

jacdet(m::AffineMap, x) = _jacdet(m, x, unsafe_matrix(m))

# If y = a*x+b, then x = inv(a)*(y-b) = inv(a)*y - inv(A)*b
inv(m::AffineMap) = AffineMap(inv(m.A), -inv(m.A)*m.b)
inverse(m::AffineMap, x) = m.A \ (x-m.b)

function leftinverse(m::AffineMap)
    M, N = size(m)
    M < N && error("No left inverse exists for $(m)")
    pA = pinv(m.A)
    AffineMap(pA, -pA*m.b)
end
function rightinverse(m::AffineMap)
    M, N = size(m)
    M > N && error("No right inverse exists for $(m)")
    pA = pinv(m.A)
    AffineMap(pA, -pA*m.b)
end
function leftinverse(m::AffineMap, x)
    M, N = size(m)
    @assert M >= N
    m.A \ (x-m.b)
end
function rightinverse(m::AffineMap, x)
    M, N = size(m)
    @assert M <= N
    m.A \ (x-m.b)
end


"An affine map for any combination of types of `A` and `b`."
struct GenericAffineMap{T,AA,B} <: AffineMap{T}
    A   ::  AA
    b   ::  B
end

GenericAffineMap(A, b) = GenericAffineMap{promote_type(eltype(A),eltype(b))}(A, b)
GenericAffineMap(A::AbstractArray{S}, b::AbstractVector{T}) where {S,T} =
    GenericAffineMap{Vector{promote_type(S,T)}}(A, b)
GenericAffineMap(A::S, b::AbstractVector{T}) where {S<:Number,T} =
    GenericAffineMap{Vector{promote_type(S,T)}}(A, b)
GenericAffineMap(A::S, b::SVector{N,T}) where {S<:Number,N,T} =
    GenericAffineMap{SVector{N,promote_type(S,T)}}(A, b)


# Fallback routine for generic A and b, special cases follow
GenericAffineMap{T}(A, b) where {T} = GenericAffineMap{T,typeof(A),typeof(b)}(A, b)

GenericAffineMap{T}(A::Number, b) where {T} = GenericAffineMap{T,eltype(T),typeof(b)}(A, b)
GenericAffineMap{T}(A::Number, b::AbstractVector) where {N,S,T <: SVector{N,S}} =
    GenericAffineMap{T,S,SVector{N,S}}(convert(S,A), SVector{N,S}(b))
# Promote element types of abstract arrays
GenericAffineMap{T}(A::AbstractMatrix, b::AbstractVector) where {S,T<:AbstractVector{S}} =
    GenericAffineMap{T}(convert(AbstractMatrix{eltype(T)},A), convert(AbstractVector{eltype(T)}, b))
GenericAffineMap{T}(A::AbstractMatrix{S}, b::AbstractVector{S}) where {S,T<:AbstractVector{S}} =
    GenericAffineMap{T,typeof(A),typeof(b)}(A, b)


convert(::Type{Map{T}}, m::GenericAffineMap{T}) where {T} = m
convert(::Type{Map{T}}, m::GenericAffineMap{S}) where {S,T} = GenericAffineMap{T}(m.A, m.b)
convert(::Type{GenericAffineMap{T}}, m::GenericAffineMap) where {T} = convert(Map{T}, m)

matrix(m::GenericAffineMap{T,<:Number}) where {T<:AbstractVector} = tomatrix(T, m.A)




"An affine map with scalar representation."
struct ScalarAffineMap{T} <: AffineMap{T}
    A   ::  T
    b   ::  T
end

ScalarAffineMap(A, b) = ScalarAffineMap(promote(A, b)...)



"An affine map with array and vector representation."
struct VectorAffineMap{T} <: AffineMap{Vector{T}}
    A   ::  Array{T}
    b   ::  Vector{T}
end

VectorAffineMap(A::AbstractArray{T}, b::AbstractVector{T}) where {T} =
    VectorAffineMap{T}(A, b)
function VectorAffineMap(A::AbstractArray{S}, b::AbstractVector{T}) where {S,T}
    U = promote_type(S,T)
    VectorAffineMap(convert(AbstractArray{U}, A), convert(AbstractVector{U}, b))
end



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
StaticAffineMap{T}(A::AbstractMatrix{T}, b::SVector{M,T}) where {T,M} =
    StaticAffineMap{T,M,M}(A, b)
StaticAffineMap{T}(A::SMatrix{M,N,T}, b::AbstractVector) where {T,N,M} =
    StaticAffineMap{T,N,M}(A, b)
StaticAffineMap{T}(A::SMatrix{M,N,T}, b::SVector{M,T}) where {T,N,M} =
    StaticAffineMap{T,N,M}(A, b)
StaticAffineMap{T,N}(A::AbstractMatrix, b::AbstractVector) where {T,N} =
    StaticAffineMap{T,N,N}(A, b)
StaticAffineMap{T,N}(A::SMatrix{M,N}, b::AbstractVector) where {T,N,M} =
    StaticAffineMap{T,N,M}(A, b)

# - finally invoke the constructor (and implicitly convert the data if necessary)
StaticAffineMap{T,N,M}(A::AbstractMatrix, b::AbstractVector) where {T,N,M} =
    StaticAffineMap{T,N,M,M*N}(A, b)

convert(::Type{Map{SVector{N,T}}}, m::VectorAffineMap) where {N,T} = StaticAffineMap{T,N}(m.A, m.b)



########################
# Some useful functions
########################

@deprecate linear_map(a, b) AffineMap(a, b)

"Map the interval [a,b] to the interval [c,d]."
interval_map(a, b, c, d) = AffineMap((d-c)/(b-a), c - a*(d-c)/(b-a))


## Simple scaling maps

"Scale all variables by a."
scaling_map(a) = LinearMap(a)

"Scale the variables by a and b."
scaling_map(a, b) = LinearMap(SMatrix{2,2}(a, 0, 0, b))

"Scale the variables by a, b and c."
scaling_map(a, b, c) = LinearMap(SMatrix{3,3}(a, 0, 0,  0, b, 0,  0, 0,c))

"Scale the variables by a, b, c and d."
scaling_map(a, b, c, d) = LinearMap(SMatrix{4,4}(a,0,0,0, 0,b,0,0, 0,0,c,0, 0,0,0,d))




##############
# Arithmetic
##############

mapcompose(m1::AbstractAffineMap, m2::AbstractAffineMap) = affine_composition(m1, m2)

"""
Compute the affine map that represents map2 after map1, that is:
`y = a2*(a1*x+b1)+b2 = a2*a1*x + a2*b1 + b2`.
"""
affine_composition(map1::AbstractAffineMap, map2::AbstractAffineMap) =
    AffineMap(matrix(map2) * matrix(map1), matrix(map2)*vector(map1) + vector(map2))

affine_composition(map1::AffineMap, map2::AffineMap) =
    AffineMap(unsafe_matrix(map2) * unsafe_matrix(map1), unsafe_matrix(map2)*unsafe_vector(map1) + unsafe_vector(map2))

affine_composition(map1::LinearMap, map2::LinearMap) =
    LinearMap(unsafe_matrix(map2) * unsafe_matrix(map1))

affine_composition(map1::LinearMap, map2::AffineMap) =
    AffineMap(unsafe_matrix(map2) * unsafe_matrix(map1), unsafe_vector(map2))

affine_composition(map1::AffineMap, map2::LinearMap) =
    AffineMap(unsafe_matrix(map2) * unsafe_matrix(map1), unsafe_matrix(map2)*unsafe_vector(map1))

affine_composition(map1::Translation, map2::Translation) =
    Translation(unsafe_vector(map2) + unsafe_vector(map1))
