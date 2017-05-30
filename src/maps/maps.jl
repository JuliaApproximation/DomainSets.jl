# maps.jl

"""
A map is any transformation of the form y = f(x).
"""
abstract type AbstractMap end

inverse_map(map::AbstractMap, y) = forward_map(inv(map), y)

(*)(map::AbstractMap, x) = forward_map(map, x)

(\)(map::AbstractMap, y) = inverse_map(map, y)

is_linear(map::AbstractMap) = false

isreal(map::AbstractMap) = true

linearize(map::AbstractMap, x) = (jacobian(map, x), translation_vector(map, x))

"""
Return the matrix and vector of a linear map, with elements of the given type
(which defaults to eltype, if applicable).
"""
function matrix_vector(map::AbstractMap, T = eltype(map))
    is_linear(map) || throw(ExceptionError())
    N = ndims(map)
    I = eye(SMatrix{N,N,T})
    B = map * zeros(SVector{N,T})
    mA = zeros(T,N,N)
    for i in 1:N
        v = I[:,i]
        mA[:,i] = map * v
    end
    A = SMatrix{N,N}(mA)
    A,B
end

is_compatible(m1::AbstractMap, m2::AbstractMap) = m1==m2
