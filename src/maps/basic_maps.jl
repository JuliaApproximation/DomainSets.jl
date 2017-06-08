# basic_maps.jl

"""
The identity map.
"""
struct IdentityMap{T,S} <: AbstractMap{T,S}
end

(m::IdentityMap)(x) = applymap(m, x)

applymap(map::IdentityMap, x) = x

apply_inverse(map::IdentityMap, y) = y

inv(::IdentityMap{T,S}) where {T,S} = IdentityMap{S,T}()

jacobian{N,T}(map::IdentityMap, x::SVector{N,T}) = eye(SMatrix{N,N,T})

jacobian{T}(map::IdentityMap, x::Vector{T}) = eye(T, length(x), length(x))

is_linear(map::IdentityMap) = true

translation_vector{N,T}(map::IdentityMap, x::SVector{N,T}) = @SVector zeros(T,N)

translation_vector{T}(map::IdentityMap, x::Vector{T}) = zeros(T,length(x))

# dest_type{T}(map::IdentityMap, ::Type{T}) = T




########################
# Special maps
########################

translation{N,T}(x::SVector{N,T}) = AffineMap(eye(SMatrix{N,N,T}), x)


# General rotational maps

"""
A 2D rotation matrix corresponding to the angle `Θ`.
"""
rotationmatrix(Θ) = SMatrix{2,2}([cos(Θ) -sin(Θ); sin(Θ) cos(Θ)])

"""
Rotation about X-axis (phi), Y-axis (theta) and Z-axis (psi).
"""
rotationmatrix(ϕ, Θ, ψ) =
    SMatrix{3,3}([cos(Θ)*cos(ψ) cos(ϕ)*sin(ψ)+sin(ϕ)*sin(Θ)*cos(ψ) sin(ϕ)*sin(ψ)-cos(ϕ)*sin(Θ)*cos(ψ);
    -cos(Θ)*sin(ψ) cos(ϕ)*cos(ψ)-sin(ϕ)*sin(Θ)*sin(ψ) sin(ϕ)*cos(ψ)+cos(ϕ)*sin(Θ)*sin(ψ);
    sin(Θ) -sin(ϕ)*cos(Θ) cos(ϕ)*cos(Θ)])

rotation(Θ) = AffineMap(rotationmatrix(Θ))

rotation(ϕ, Θ, ψ) = AffineMap(rotationmatrix(ϕ,Θ,ψ))
