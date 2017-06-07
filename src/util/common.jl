# common.jl


#######################
# Composite structures
#######################

"""
Some types have composite structure, e.g. product domains, a union of domains.
These types contain a list of domains.

It is often undesirable to use `getindex` to access the elements of the composite
type. For this reason we introduce the `elements` functions. Composite types
can implement `elements` and provide a generic way to access their components.

`elements(t)`: returns the elements making up the composite type `t`

`element(t, i)`: return the `i`-th element of the composite type `t`

`nb_elements(t)`: return the number of elements of the composite type `t`
"""
# Forward declaration
elements() = nothing

"""
Return the i-th element of a composite structure.

See also: `elements`.
"""
# By default, we index elements(o)
element(t, i) = elements(t)[i]

"""
Return the number of elements of a composite structure.

See also: `elements`.
"""
# By default, we return length(elements(t))
nb_elements(t) = length(elements(t))


###############
# Type factory
###############

"""
A `TypeFactory{T}` is a convenience type to simplify construction of a type.

Having `t = TypeFactory{T}` overrides `getindex` such that `t[a]` invokes `T(a)`.

For example:
```
v = TypeFactory{SVector}
v[0.1,0.2]
```
makes an `SVector{2,Float64}`.
"""
struct TypeFactory{T}
end

Base.getindex(v::TypeFactory{T}, args...) where {T} = T(args...)

# const v = TypeFactory{SVector}()
