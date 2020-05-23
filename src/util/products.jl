
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


"""
Convert a vector from a cartesian format to a nested tuple according to the
given dimensions.

For example:
`convert_fromcartesian([1,2,3,4,5], Val{(2,2,1)}()) -> ([1,2],[3,4],5)`
"""
@generated function convert_fromcartesian(x::AbstractVector, ::Val{DIM}) where {DIM}
	dimsum = [0; cumsum([d for d in DIM])]
	E = Expr(:tuple, [ (dimsum[i+1]-dimsum[i] > 1 ? Expr(:call, :SVector, [:(x[$j]) for j = dimsum[i]+1:dimsum[i+1]]...) : :(x[$(dimsum[i+1])])) for i in 1:length(DIM)]...)
	return quote $(E) end
end

"The inverse function of `convert_fromcartesian`."
@generated function convert_tocartesian(x, ::Val{DIM}) where {DIM}
    dimsum = [0; cumsum([d for d in DIM])]
    E = vcat([[:(x[$i][$j]) for j in 1:DIM[i]] for i in 1:length(DIM)]...)
    quote SVector($(E...)) end
end
# An alternative is to use "reduce(vcat, x)" (see Julia issue #21672) but the
# generated function is more efficient because the compiler knows the dimensions.
