# common.jl

# Forward declarations
"""
Return all the elements of a composite structure.
"""
function elements() end

"""
Return the i-th element of a composite structure.
"""
element(s, i) = elements(s)[i]

"""
Return the number of elements of a composite structure.
"""
nb_elements(s) = length(elements(s))

# for legacy reasons - remove soon
composite_length = nb_elements
