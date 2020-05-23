# Routines having to do with computations involving domains

@deprecate +(d::Domain, x::Union{AbstractVector,Number}) d .+ x
@deprecate -(d::Domain, x::Union{AbstractVector,Number}) d .- x
@deprecate +(x::Union{Number,AbstractVector}, d::Domain) x .+ d
@deprecate -(x::Union{Number,AbstractVector}, d::Domain) x .- d

# Allow unary minus, but use broadcast for the implementation
-(d::Domain) = (-).(d)

# Allow multiplication and division by numbers, like for vectors
*(a::Number, domain::Domain) = a .* domain
*(domain::Domain, a::Number) = domain .* a
/(domain::Domain, a::Number) = domain ./ a
\(a::Number, domain::Domain) = a .\ domain
