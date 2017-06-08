# coordinates.jl
# Maps having to do with coordinate transforms.

"""
A Cartesion to Polar map. First dimension is interpreted as angle in radians, second as radial distance
A square [-1,1]x[-1,1] is mapped to the unit circle

"""
struct CartToPolarMap{T} <: AbstractMap{SVector{2,T}, SVector{2,T}}
end

(m::CartToPolarMap)(x) = applymap(m, x)

applymap(map::CartToPolarMap, x::SVector{2,T}) where {T} = SVector{2,T}(((x[2]+1)/2)*cos(pi*x[1]), ((x[2]+1)/2)*sin(pi*x[1]))

inv(map::CartToPolarMap{T}) where {T} = PolarToCartMap{T}()

is_linear(map::CartToPolarMap) = false


"""
A Polar to Cartesian map. The angle is mapped to the first dimension, radius to the second.
The unit circle is mapped to a square [-1,1]x[-1,1]
"""
struct PolarToCartMap{T} <: AbstractMap{SVector{2,T}, SVector{2,T}}
end

(m::PolarToCartMap)(x) = applymap(m, x)

applymap(map::PolarToCartMap, x::SVector{2,T}) where {T} = SVector{2,T}(sqrt(x[1]^2+x[2]^2)*2-1, atan2(x[2],x[1])/pi)

inv(map::PolarToCartMap{T}) where {T} = CartToPolarMap{T}()

is_linear(map::PolarToCartMap) = false
