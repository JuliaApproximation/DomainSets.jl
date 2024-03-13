
# Simplifications go here

# Basic maps

composedmap1(m1::IdentityMap, m2) = m2
composedmap1(m1::IdentityMap{T}, m2::Map{T}) where {T} = m2
composedmap1(m1::IdentityMap{T}, m2::Map{S}) where {S,T} = convert(Map{T}, m2)

composedmap2(m1, m2::IdentityMap) = m1
composedmap2(m1::Map{T}, m2::IdentityMap{T}) where {T} = m1
composedmap2(m1::Map{S}, m2::IdentityMap{T}) where {S,T} = convert(Map{T}, m1)

composedmap2(m1, m2::ConstantMap) = m2
composedmap2(m1::Map{T}, m2::ConstantMap) where {T} = ConstantMap{T}(mapconstant(m2))
composedmap1(m1::ConstantMap{T}, m2) where {T} = ConstantMap{T}(m2(mapconstant(m1)))

composedmap2(m1, m2::ZeroMap) = m2
composedmap2(m1::Map{T}, m2::ZeroMap{S,U}) where {S,T,U} = ZeroMap{T,U}()
composedmap1(m1::ConstantMap{T}, m2::ZeroMap{S,U}) where {S,T,U} = ZeroMap{T,U}()

multiply_map1(m1::ZeroMap, m2) = m1
multiply_map2(m1, m2::ZeroMap) = m2
multiply_map2(m1, m2::ConstantMap) = _constant_multiply_map(m1, m2)
_constant_multiply_map(m1::ConstantMap{T}, m2::ConstantMap{S}) where {S,T} =
    ConstantMap{promote_type(S,T)}(mapconstant(m1)*mapconstant(m2))
_constant_multiply_map(m1, m2::ConstantMap) = default_multiply_map(m1, m2)

multiply_map(m1::Function, m2::Function) = t -> m1(t)*m2(t)

sum_map1(m1::ZeroMap, m2) = m2
sum_map2(m1, m2::ZeroMap) = m1
sum_map2(m1, m2::ConstantMap) = _constant_sum_map(m1, m2)
_constant_sum_map(m1::ConstantMap{T}, m2::ConstantMap{S}) where {S,T} =
    ConstantMap{promote_type(S,T)}(mapconstant(m1)+mapconstant(m2))
_constant_sum_map(m1, m2::ConstantMap) = default_sum_map(m1, m2)

## Affine maps

composedmap(m1::AbstractAffineMap, m2::AbstractAffineMap) = affine_composition(m1, m2)

"""
Compute the affine map that represents map2 after map1, that is:
`y = a2*(a1*x+b1)+b2 = a2*a1*x + a2*b1 + b2`.
"""
affine_composition(map1::AbstractAffineMap, map2::AbstractAffineMap) =
    AffineMap(affinematrix(map2) * affinematrix(map1),
        affinematrix(map2)*affinevector(map1) + affinevector(map2))

affine_composition(map1::AffineMap, map2::AffineMap) =
    AffineMap(unsafe_matrix(map2) * unsafe_matrix(map1),
        unsafe_matrix(map2)*unsafe_vector(map1) + unsafe_vector(map2))

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
sum_map(map1::AbstractAffineMap, map2::AbstractAffineMap) =
    AffineMap(affinematrix(map1)+affinematrix(map2), affinevector(map1)+affinevector(map2))


isequalmap(m1::ProductMap, m2::IdentityMap) = all(map(isidentitymap, components(m1)))
isequalmap(m1::IdentityMap, m2::ProductMap) = isequalmap(m2, m1)


"""
Map the interval `[a,b]` to the interval `[c,d]`.

This function deals with infinite intervals, and the type of the
map returned may depend on the value (finiteness) of the given endpoints.
"""
interval_map(a, b, c, d) = interval_map(promote(a,b,c,d)...)

function interval_map(a::T, b::T, c::T, d::T) where {T}
    FT = float(T)
    if isfinite(a) && isfinite(b) && isfinite(c) && isfinite(d)
        bounded_interval_map(a, b, c, d)
    elseif isfinite(a) && !isfinite(b) && isfinite(c) && !isfinite(d)
        # (a,Inf) to (c,Inf)
        AffineMap(one(FT), c-a)
    elseif isfinite(a) && !isfinite(b) && !isfinite(c) && isfinite(d)
        # (a,Inf) to (Inf,d)
        AffineMap(-one(FT), d+a)
    elseif !isfinite(a) && isfinite(b) && isfinite(c) && !isfinite(d)
        # (Inf,b) to (c,Inf)
        AffineMap(-one(FT), c+b)
    elseif !isfinite(a) && isfinite(b) && !isfinite(c) && isfinite(d)
        # (Inf,b) to (Inf,d)
        AffineMap(one(FT), d-b)
    elseif !isfinite(a) && !isfinite(b) && !isfinite(c) && !isfinite(d)
        if (a < 0) && (b > 0) && (c < 0) && (d > 0)
            # (-Inf,Inf) to (-Inf,Inf)
            StaticIdentityMap{FT}()
        elseif (a < 0) && (b > 0) && (c > 0) && (d < 0)
            # (-Inf,Inf) to (Inf,-Inf)
            LinearMap(-one(FT))
        elseif (a > 0) && (b < 0) && (c < 0) && (d > 0)
            # (Inf,-Inf) to (-Inf,Inf)
            LinearMap(-one(FT))
        elseif (a > 0) && (b < 0) && (c > 0) && (d < 0)
            # (Inf,-Inf) to (Inf,-Inf)
            StaticIdentityMap{FT}()
        elseif (a > 0) && (b > 0) && (c > 0) && (d > 0)
            # (Inf,Inf) to (Inf,Inf)
            StaticIdentityMap{FT}()
        elseif (a < 0) && (b < 0) && (c < 0) && (d < 0)
            # (-Inf,-Inf) to (-Inf,-Inf)
            StaticIdentityMap{FT}()
        elseif (a < 0) && (b < 0) && (c > 0) && (d > 0)
            # (-Inf,-Inf) to (Inf,Inf)
            LinearMap(-one(FT))
        elseif (a > 0) && (b > 0) && (c < 0) && (d < 0)
            # (Inf,Inf) to (-Inf,-Inf)
            LinearMap(-one(FT))
        else
            throw(ArgumentError("Requested affine map is unbounded"))
        end
    else
        throw(ArgumentError("Requested affine map is unbounded"))
    end
end

"Like interval_map, but guaranteed to return a scalar affine map."
bounded_interval_map(a, b, c, d) = bounded_interval_map(promote(a,b,c,d)...)
bounded_interval_map(a::T, b::T, c::T, d::T) where {T} =
    AffineMap((d-c)/(b-a), c - a*(d-c)/(b-a))
