# common.jl

# Forward declarations
"""
Return all the elements of a composite structure.
"""
elements() = nothing

"""
Return the i-th element of a composite structure.
"""
element(s, i) = elements(s)[i]

"""
Return the number of elements of a composite structure.
"""
nb_elements(s) = length(elements(s))
