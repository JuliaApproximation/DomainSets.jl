# In this file we group a number of concepts relating to geometrical features.
# These routines may make assumptions about orientations and a metric that
# might not be entirely generic.

"""
Return a bounding box of the given domain.

A bounding box is an interval, a hyperrectangle or the full space. It is such that
each point in the domain also lies in the bounding box.
"""
boundingbox(d) = fullspace(d)

"Return the boundary of the given domain as a domain."
function boundary end
const ∂ = boundary

"""
Return the normal of the domain at the point `x`.

It is assumed that `x` is a point on the boundary of the domain.
"""
function normal end

"""
Return the tangents of the domain at the point `x`. The tangents form a
basis for the tangent plane, perpendicular to the normal direction at `x`.
"""
function tangents end

"""
    distance_to(d, x)

Return the distance from the point `x` to the domain `d`.
"""
function distance_to end

"Return the volume of the domain."
function volume end
