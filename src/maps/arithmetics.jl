
# Simplifications go here

compose_map1(m1::IdentityMap, m2) = m2
compose_map2(m1, m2::IdentityMap) = m1

compose_map(m1::AbstractAffineMap, m2::AbstractAffineMap) = affine_composition(m1, m2)

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

affine_composition(map1::Translation, map2::LinearMap) =
    AffineMap(unsafe_matrix(map2), unsafe_matrix(map2)*unsafe_vector(map1))

affine_composition(map1::LinearMap, map2::Translation) =
    AffineMap(unsafe_matrix(map1), unsafe_vector(map2))

# The sum of two affine maps is again an affine map
mapsum(map1::AbstractAffineMap, map2::AbstractAffineMap) =
    AffineMap(matrix(map1)+matrix(map2), vector(map1)+vector(map2))


==(m1::ProductMap, m2::IdentityMap) = all(map(isidentity, components(m1)))
==(m1::IdentityMap, m2::ProductMap) = m2 == m1
