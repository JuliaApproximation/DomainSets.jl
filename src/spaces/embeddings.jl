# embeddings.jl


#######################################
# Isomorphism between geometric spaces
#######################################

"""
True if the two geometric spaces are isomorphic.

If they are, then a point in one space can be converted into
a point in the other space and vice-versa, using the `convert_space` function.
"""
# By default spaces are isomorphic only if they are identical
isomorphic(s1::GeometricSpace, s2::GeometricSpace) = (s1 == s2)

"""
Convert the point `x` from the space `s1` to a point in the space `s2`.
This is only possible if the spaces are isomorphic.
"""
function convert_space(x, s1::GeometricSpace, s2::GeometricSpace)
    @assert isomorphic(s1, s2)
    @assert x ∈ s1
    if s1 == s2
        # Don't do anything if the spaces are identical
        x
    else
        # Concrete subtypes have to implement unsafe_convert_space
        unsafe_convert_space(x, s1, s2)
    end
end

# We can provide a default based on the Julia conversion system
unsafe_convert_space(x, s1::GeometricSpace, s2::GeometricSpace) = convert(eltype(s2), x)


#################################
# Embedding of geometric spaces
#################################

"""
True if the first geometric space is embedded into the second space.

If it is, then a point in `s1` can be promoted to a point in `s2`, using the
`promote_space` function. A point in `s2` can be restricted to `s1`, if it lies
in the image of the canonical embedding, using `restrict_space`. All points in
`s2` can be projected onto points in `s1` using `project`.
"""
# By default, isomorphic spaces are embedded
embedded(s1::GeometricSpace, s2::GeometricSpace) = isomorphic(s1, s2)

"""
Promote the point `x` from the space `s1` to a point in the space `s2`.
This is only possible if the `s1` is embedded into `s2`.
"""
function promote_space(x, s1::GeometricSpace, s2::GeometricSpace)
    @assert embedded(s1, s2)
    @assert x ∈ s1
    if s1 == s2
        # Don't do anything if the spaces are identical
        x
    else
        # Concrete subtypes have to implement unsafe_promote_space
        unsafe_promote_space(x, s1, s2)
    end
end

# We can provide a default based on the Julia conversion system
unsafe_promote_space(x, s1::GeometricSpace, s2::GeometricSpace) = convert(eltype(s2), x)
