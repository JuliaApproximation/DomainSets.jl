# Rotation around the origin
rotate(d::EuclideanDomain{2}, θ) = rotation_map(θ).(d)

rotate(d::EuclideanDomain{3}, phi, theta, psi) = rotation_map(phi,theta,psi).(d)
# Rotation around a fixed center.
rotate(d::EuclideanDomain{2}, θ, center::SVector{T}) where {T} = (Translation(center) ∘ rotation_map(θ) ∘ Translation(-center)).(d)

rotate(d::EuclideanDomain{3}, phi, theta, psi, center::SVector{T}) where {T} = (Translation(center) ∘ rotation_map(phi,theta,psi) ∘ Translation(-center)).(d)


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
