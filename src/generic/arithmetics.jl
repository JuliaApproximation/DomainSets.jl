# Routines having to do with computations involving domains

@deprecate +(d::Domain, x::Union{AbstractVector,Number}) d .+ x
@deprecate -(d::Domain, x::Union{AbstractVector,Number}) d .- x
@deprecate +(x::Union{Number,AbstractVector}, d::Domain) x .+ d
@deprecate -(x::Union{Number,AbstractVector}, d::Domain) x .- d
@deprecate -(d::Domain) (-).(d)
@deprecate *(a::Number, domain::Domain) a .* domain
@deprecate *(domain::Domain, a::Number) domain .* a
@deprecate /(domain::Domain, a::Number) domain ./ a
@deprecate \(a::Number, domain::Domain) a .\ domain
