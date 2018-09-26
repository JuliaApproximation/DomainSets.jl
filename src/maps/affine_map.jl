# affine_map.jl

"""
An affine map has the general form `y = a*x + b`, with types for `a`, `b`, `x`
and `y` such that the expression is valid.

We use matrix and vector to denote `a` and `b` respectively.
"""
abstract type AbstractAffineMap{S,T} <: AbstractMap{S,T}
end

jacobian(m::AbstractAffineMap, x) = matrix(m)

function jacobian(m::AbstractAffineMap{S,T}) where {S,T}
    A = matrix(m) * one(jac_type(S,T))
    ConstantMap{S,typeof(A)}(A)
end

islinear(map::AbstractMap) = false
islinear(map::AbstractAffineMap) = true

update_eltype(map::AbstractAffineMap{S,T}, ::Type{T}) where {S,T} = map

update_eltype(map::AbstractAffineMap{S,T}, ::Type{U}) where {S,T,U} =
    map_update_eltype(map, U)

"""
A `LinearMap` is an affine map that represents `y = a*x`, where `a` can have any
type such that `a*x` maps type `S` to type `T.`
"""
struct LinearMap{S,T,A} <: AbstractAffineMap{S,T}
    a   ::  A
end

LinearMap{T}(a) where {T} = LinearMap{T,T}(a)

LinearMap{S,T}(a) where {S,T} = LinearMap{S,T,typeof(a)}(a)

LinearMap(a::SMatrix{M,N,T}) where {M,N,T} = LinearMap{SVector{N,T},SVector{M,T},typeof(a)}(a)

LinearMap(a::Matrix{T}) where {T} = LinearMap{Vector{T},Vector{T},typeof(a)}(a)

LinearMap(a::T) where {T <: Number} = LinearMap{T,T,T}(a)

matrix(m::LinearMap) = m.a

vector(m::LinearMap) = zero(codomaintype(m))

function map_update_eltype(m::LinearMap, T)
    a = map(T, matrix(m))
    LinearMap(a)
end

(m::LinearMap)(x) = applymap(m, x)

applymap(m::LinearMap, x) = matrix(m) * x

inv(m::LinearMap{S,T}) where {S,T} = LinearMap{T,S}(inv(matrix(m)))

# Because StaticArrays does not currently support `pinv` we include a workaround:
pinv(m::SMatrix{M,N}) where {M,N}  = SMatrix{N,M}(pinv(convert(Array,m)))

left_inverse(m::LinearMap{S,T}) where {S,T} =  LinearMap{T,S}(pinv(matrix(m)))
right_inverse(m::LinearMap{S,T}) where {S,T} = LinearMap{T,S}(pinv(matrix(m)))


"""
Translation represents `y = x + v`, where `v` is a vector in the same space as
`x` and `y`.
"""
struct Translation{T} <: AbstractAffineMap{T,T}
    vector  ::  T
end

translation_map(vector::T) where T = Translation{T}(vector)

matrix(m::Translation{T}) where T = one(jac_type(T,T))

vector(m::Translation) = m.vector

function map_update_eltype(m::Translation, T)
    b = map(T, vector(m))
    Translation(b)
end

(m::Translation)(x) = applymap(m, x)

applymap(m::Translation, x) = x + vector(m)

apply_inverse(m::Translation, y) = y - vector(m)

inv(m::Translation) = Translation(-vector(m))


"""
`AffineMap` represents `y = a*x + b`, i.e. it combines a `LinearMap` and a
`Translation`.
"""
struct AffineMap{S,T,A} <: AbstractAffineMap{S,T}
    a   ::  A
    b   ::  T
end

AffineMap{S,T}(a::A, b) where {S,T,A} = AffineMap{S,T,A}(a, b)

AffineMap(a::T, b::T) where {T} = AffineMap{T,T,T}(a, b)

AffineMap(a::SMatrix{M,N,T}, b::SVector{M,T}) where {M,N,T} =
    AffineMap{SVector{N,T},SVector{M,T},typeof(a)}(a, b)

AffineMap(a::Number, b::SVector{N,T}) where {N,T} =
     AffineMap(SMatrix{N,N,T}(1.0I), b)

matrix(m::AffineMap) = m.a

vector(m::AffineMap) = m.b

function map_update_eltype(m::AffineMap, T)
    a = map(T, matrix(m))
    b = map(T, vector(m))
    AffineMap(a, b)
end

(m::AffineMap)(x) = applymap(m, x)

applymap(m::AffineMap, x) = m.a * x .+ m.b

# If y = a*x+b, then x = inv(a)*(y-b).
inv(m::AffineMap{S,T}) where {S,T} = AffineMap{T,S}(inv(m.a), -inv(m.a)*m.b)

# Todo: implement `full` for affine maps, which would result in `a` always being
# a dense matrix.


########################
# Some useful functions
########################

"Make the linear map y = a*x + b."
linear_map(a, b) = AffineMap(a, b)

"Map the interval [a,b] to the interval [c,d]."
interval_map(a, b, c, d) = linear_map((d-c)/(b-a), c - a*(d-c)/(b-a))


## Simple scaling maps

"Scale all variables by a."
scaling_map(a) = LinearMap(a)

"Scale the variables by a and b."
scaling_map(a, b) = LinearMap(SMatrix{2,2}(a, 0, 0, b))

"Scale the variables by a, b and c."
scaling_map(a, b, c) = LinearMap(SMatrix{3,3}(a, 0, 0,  0, b, 0,  0, 0,c))

"Scale the variables by a, b, c and d."
scaling_map(a, b, c, d) = LinearMap(SMatrix{4,4}(a,0,0,0, 0,b,0,0, 0,0,c,0, 0,0,0,d))


#############################
# Rotations around the origin
#############################

# Rotation in positive (counterclockwise) direction
# (Note: the SMatrix constructor expects the arguments column-first)
rotationmatrix(theta) = SMatrix{2,2}(cos(theta), sin(theta), -sin(theta), cos(theta))

# Rotation about X-axis (phi), Y-axis (theta) and Z-axis (psi)
# As above, the matrix is given column-by-column
rotationmatrix(phi,theta,psi) =
    SMatrix{3,3}(cos(theta)*cos(psi), cos(theta)*sin(psi), -sin(theta),
        -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi), cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi), sin(phi)*cos(theta),
        sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi), -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi), cos(phi)*cos(theta))

rotation_map(theta) = LinearMap(rotationmatrix(theta))

rotation_map(phi, theta, psi) = LinearMap(rotationmatrix(phi,theta,psi))


##############
# Arithmetic
##############

∘(m2::AbstractAffineMap, m1::AbstractAffineMap) = affine_composition(m2, m1)

"""
Compute the affine map that represents map2 after map1, that is:
y = a2*(a1*x+b1)+b2 = a2*a1*x + a2*b1 + b2.
"""
affine_composition(map2::AbstractAffineMap, map1::AbstractAffineMap) =
    AffineMap(matrix(map2) * matrix(map1), matrix(map2)*vector(map1) + vector(map2))

affine_composition(map2::LinearMap{U,T}, map1::LinearMap{T,S}) where {S,T,U} =
    LinearMap{S,U}(matrix(map2) * matrix(map1))

affine_composition(map2::Translation{T}, map1::Translation{T}) where {T} =
    Translation{T}(translation_vector(map2) + translation_vector(map1))
