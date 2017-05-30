# affine_map.jl

"""
An affine map has the form `y = a*x + b`.

The fields a and b can be anything that can multiply a vector (a) and be added
to a vector (b). In higher dimensions, the intended use is for a and b to be
static arrays. In that case, the application of the map does not allocate memory.

The affinemap also stores its inverse `x = c*y + d`, with `c = inv(a)` and
`d = - inv(a)*b`. It is enforced that `a` and `c` have the same type, as do
`b` and `d`.

This restricts the possibilities for `a` and `b` somewhat. For example, `b` can
not be a number if `a` is a matrix, since in that case `-inv(a)*b` would be a
vector. An exception is when `b` is exactly `0`.

In general, `a` is like a matrix and `b` is like a vector. However, there are
other valid cases:
* `a` is a scalar and `b` is a vector
* `a` is a matrix and `b` is exactly 0
"""
struct AffineMap{TA,TB} <: AbstractMap
    # The fields a and b define the forward map y = a*x+b
    a   ::  TA
    b   ::  TB
    # The fields c and d define the inverse map x = c*y+d
    c   ::  TA
    d   ::  TB
end

# If only one argument is given, we assume that b is 0.
AffineMap(a) = AffineMap(a, 0)

# If only two arguments are given we compute the inverse.
function AffineMap(a, b)
    c, d = affine_inv(a, b)
    AffineMap(a, b, c, d)
end

# This outer constructor will be called if the types of a and c or b and d do
# not match. We promote and call the inner constructor.
function AffineMap(a, b, c, d)
    a, c = promote(a, c)
    b, d = promote(b, d)
    AffineMap(a, b, c, d)
end

## The logic for the inverse of a*x+b follows.

"Compute c and d such that x = c*y + d is the inverse of y = a*x + b"
affine_inv(a, b) = (inv(a), -a\b)

# If b==0 then we can get by with b = 0 as well.
function affine_inv(a::AbstractMatrix, b::Number)
    @assert b == 0
    (inv(a), 0)
end


(m::AffineMap)(x) = forward_map(m, x)

eltype(map::AffineMap) = promote_type(eltype(map.a), eltype(map.b))

forward_map(map::AffineMap, x) = map.a * x + map.b

inverse_map(map::AffineMap, y) = map.c * y + map.d

inv(map::AffineMap) = AffineMap(map.c, map.d, map.a, map.b)

isreal(map::AffineMap) = isreal(map.a) && isreal(map.b)


jacobian(map::AffineMap, x) = map.a

is_linear(map::AffineMap) = true

translation_vector(map::AffineMap, x) = map.b - map.a*x


# Some useful functions
"Make the linear map y = a*x + b."
linear_map(a, b) = AffineMap(a, b)

"Map the interval [a,b] to the interval [c,d]."
interval_map(a, b, c, d) = linear_map((d-c)/(b-a), c - a*(d-c)/(b-a))

## Simple scailng maps

"Scale all variables by a."
scaling_map(a) = AffineMap(a)

"Scale the variables by a and b."
scaling_map(a, b) = AffineMap(SMatrix{2,2}(a,0,0,b))

"Scale the variables by a, b and c."
scaling_map(a, b, c) = AffineMap(SMatrix{3,3}(a,0,0, 0,b,0, 0,0,c))

# 4x4 StaticArrays don't seem to work as well as 1-3d, the code below errors
# "Scale the variables by a, b, c and d."
# scaling_map(a, b, c, d) = AffineMap(SMatrix{4,4}(a,0,0,0, 0,b,0,0, 0,0,c,0, 0,0,0,d))

"""
Compute the affine map that represents map2*map1, that is:
y = a2*(a1*x+b1)+b2 = a2*a1*x + a2*b1 + b2.
"""
affine_composition(map1::AffineMap, map2::AffineMap) = affine_composition(map1.a, map1.b, map2.a, map2.b)

# We have to compute a matrix a2*a1 and a vector a2*b1+b2.
# We have to be careful to treat the cases where a and/or b are scalars properly,
# since a*x+b sometimes relies on broadcasting.
# It turns out the only problematic case is a*b when a is a matrix and b is a
# scalar. In that case a*b is again a matrix, but we want it to be a vector.
# But given the restrictions on types, we can only have this case when b = 0.
affine_composition(a1, b1, a2, b2) = AffineMap(a2*a1, affine_composition_vector(a2, b1, b2))

# About the composition vector:
# - The general expression is fine in most cases
affine_composition_vector(a2, b1, b2) = a2*b1 + b2

# - be careful when a2 is a matrix and b1 a scalar
function affine_composition_vector(a2::AbstractMatrix, b1::Number, b2)
    @assert b1 == 0
    b2
end


########################
# Arithmetic
########################

(*)(map1::AffineMap, map2::AffineMap) = affine_composition(map2, map1)

(*)(map1::AbstractMap, map2::AffineMap) = composite_map(map1, map2, is_linear(map1))

function composite_map(map1::AbstractMap, map2::AffineMap, islinear::Bool)
    if islinear
        CompositeMap((map2,map1))
    else
        A,B = matrix_vector(map1)
        AffineMap(A,B) * map2
    end
end

(*)(a::Number, m::AbstractMap) = scaling_map(a) * m

(+)(m::AffineMap, x) = AffineMap(m.a, m.b+x)

(+)(m1::AffineMap, m2::AffineMap) = AffineMap(m1.a+m2.a, m1.b+m2.b)
