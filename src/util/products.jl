# products.jl

# Functions having to do with the creation of cartesian products.

"""
Create a cartesian product of the supplied arguments.

The function cartesianproduct applies some simplifications and does not necessarily
return a Product type.

A `cartesianproduct(a)` with just a single element returns `a`.

For integer `n`, `cartesianproduct(a, n)` becomes `cartesianproduct(a, a, ..., a)`.
A type-safe variant is `cartesianproduct(a, Val{N})`.
"""
cartesianproduct() = nothing

# Don't create a cartesian product of just one element
cartesianproduct(a::Tuple) = cartesianproduct(a...)
cartesianproduct(a) = a

# Create a cartesian product with n times the same element
cartesianproduct(a, n::Int) = (n==1) ? a : cartesianproduct(ntuple(t->a, n)...)

cartesianproduct(a, ::Type{Val{N}}) where {N} = cartesianproduct(ntuple(t->a, Val(N))...)
