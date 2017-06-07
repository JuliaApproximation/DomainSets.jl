# tensorproducts.jl

# Functions having to do with the creation of tensor products.

"""
Create a tensor product of the supplied arguments.

The function tensorproduct applies some simplifications and does not necessarily
return a Product type.

A `tensorproduct(a)` with just a single element returns `a`.

For integer `n`, `tensorproduct(a, n)` becomes `tensorproduct(a, a, ..., a)`.
A type-safe variant is `tensorproduct(a, Val{N})`.
"""
tensorproduct() = nothing

# Don't create a tensor product of just one element
tensorproduct(a::Tuple) = tensorproduct(a...)
tensorproduct(a) = a

# Create a tensor product with n times the same element
tensorproduct(a, n::Int) = tensorproduct(ntuple(t->a, n)...)

tensorproduct(a, ::Type{Val{N}}) where {N} = tensorproduct(ntuple(t->a, Val{N})...)

# Use \otimes as notation for tensor product.
⊗ = tensorproduct
