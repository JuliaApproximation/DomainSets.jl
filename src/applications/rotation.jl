# Rotation around the origin
rotate(d::EuclideanDomain{2}, θ) = rotation_map(θ).(d)

rotate(d::EuclideanDomain{3}, phi, theta, psi) = rotation_map(phi,theta,psi).(d)
# Rotation around a fixed center.
rotate(d::EuclideanDomain{2}, θ, center::SVector{T}) where {T} = (Translation(center) ∘ rotation_map(θ) ∘ Translation(-center)).(d)

rotate(d::EuclideanDomain{3}, phi, theta, psi, center::SVector{T}) where {T} = (Translation(center) ∘ rotation_map(phi,theta,psi) ∘ Translation(-center)).(d)

# Maps having to do with coordinate transforms.

"""
A Cartesion to Polar map. First dimension is interpreted as radial distance,
second as an angle.
The unit circle is mapped to the square `[-1,1]x[-1,1]`.
"""
struct CartToPolarMap{T} <: Map{SVector{2,T}}
end

CartToPolarMap() = CartToPolarMap{Float64}()

mapsize(m::CartToPolarMap) = (2,2)

applymap(map::CartToPolarMap{T}, x) where {T} =
    SVector{2,T}(sqrt(x[1]^2+x[2]^2)*2-1, atan(x[2],x[1])/pi)

function jacobian(m::CartToPolarMap{T}, x) where {T}
    d = sqrt(x[1]^2+x[2]^2)
    r = 1/(1+(x[2]/x[1])^2)
    SMatrix{2,2,T}(2*x[1]/d, 1/pi*r*(-x[2]/x[1]^2), 2*x[2]/d, 1/pi*r*1/x[1])
end

inverse(m::CartToPolarMap{T}) where {T} = PolarToCartMap{T}()
inverse(m::CartToPolarMap, x) = inverse(m)(x)

isreal(m::CartToPolarMap) = true

convert(::Type{Map{SVector{2,T}}}, ::CartToPolarMap) where {T} = CartToPolarMap{T}()

isequalmap(m1::CartToPolarMap, m2::CartToPolarMap) = true


"""
A Polar to Cartesian map. The angle is mapped to the second dimension,
radius to the first.
The square `[-1,1]x[-1,1]` is mapped to the unit circle.
"""
struct PolarToCartMap{T} <: Map{SVector{2,T}}
end

PolarToCartMap() = PolarToCartMap{Float64}()

mapsize(m::PolarToCartMap) = (2,2)

applymap(map::PolarToCartMap{T}, x) where {T} = SVector{2,T}((x[1]+1)/2*cos(pi*x[2]), (x[1]+1)/2*sin(pi*x[2]))

jacobian(m::PolarToCartMap{T}, x) where {T} =
    SMatrix{2,2,T}(cos(pi*x[2])/2, sin(pi*x[2])/2, -pi*(x[1]+1)/2*sin(pi*x[2]), pi*(x[1]+1)/2*cos(pi*x[2]))

inverse(m::PolarToCartMap{T}) where {T} = CartToPolarMap{T}()
inverse(m::PolarToCartMap, x) = inverse(m)(x)

isreal(m::PolarToCartMap) = true

convert(::Type{Map{SVector{2,T}}}, ::PolarToCartMap) where {T} = PolarToCartMap{T}()

isequalmap(m1::PolarToCartMap, m2::PolarToCartMap) = true


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
