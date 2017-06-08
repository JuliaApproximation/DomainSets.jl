# maps.jl

"""
A map is any transformation of the form `y = f(x)`, where `x` has type `S` and
`y` has type `T`.
"""
abstract type AbstractMap{T,S}
end

domaintype(map::AbstractMap{T,S}) where {T,S} = S
rangetype(map::AbstractMap{T,S}) where {T,S} = T

(*)(map::AbstractMap, x) = applymap(map, x)

(\)(map::AbstractMap, y) = inv(map) * y

"""
`return_type(map, U)` is a generic function that computes the return type when
the given map is applied to a variable of type `U`.

For any `AbstractMap{T,S}`, we have that `return_type(map, S) = T`. We also return
`T` if `U` can be promoted to `T` using Julia's promotion system.
"""
# - If the map maps S to T and x has type S, then the result is T
return_type(map::AbstractMap{T,S}, ::Type{S}) where {T,S} = T
# - If x has a type different from S, perhaps it can be promoted to S?
return_type(map::AbstractMap{T,S}, ::Type{U}) where {T,S,U} = _return_type(map, T, promote_to(U,S))
# -- Yes, it can, let's return T.
_return_type(map, ::Type{T}, ::Type{True}) where {T} = T
# -- No, it can not. We return Any.
_return_type(map, ::Type{T}, ::Type{False}) where {T} = Any


# is_linear(map::AbstractMap) = false
#
# isreal(map::AbstractMap) = true
#
# linearize(map::AbstractMap, x) = (jacobian(map, x), translation_vector(map, x))
#
# """
# Return the matrix and vector of a linear map, with elements of the given type
# (which defaults to eltype, if applicable).
# """
# function matrix_vector(map::AbstractMap, T = eltype(map))
#     is_linear(map) || throw(ExceptionError())
#     N = ndims(map)
#     I = eye(SMatrix{N,N,T})
#     B = map * zeros(SVector{N,T})
#     mA = zeros(T,N,N)
#     for i in 1:N
#         v = I[:,i]
#         mA[:,i] = map * v
#     end
#     A = SMatrix{N,N}(mA)
#     A,B
# end
#
# is_compatible(m1::AbstractMap, m2::AbstractMap) = m1==m2
