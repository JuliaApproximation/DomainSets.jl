# tensorproducts

"""
Create a tensor product of the supplied arguments.
"""
tensorproduct() = nothing

# Use \otimes as notation for tensor product.
⊗ = tensorproduct

# Flatten a sequence of elements that may be recursively composite
# For example: a ProductDomain of ProductDomains will yield a list of each of the
# individual domains, like the leafs of a tree structure.
function flatten{T}(::Type{T}, elements::Array, BaseType = Any)
    flattened = BaseType[]
    for element in elements
        append_flattened!(T, flattened, element)
    end
    flattened
end

flatten{T}(::Type{T}, elements...) = tuple(flatten(T, [el for el in elements])...)

function append_flattened!{T}(::Type{T}, flattened::Vector, element::T)
    for el in elements(element)
        append_flattened!(T, flattened, el)
    end
end

function append_flattened!{T}(::Type{T}, flattened::Vector, element)
    append!(flattened, [element])
end
