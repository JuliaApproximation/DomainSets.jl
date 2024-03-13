
# some renames, to be deprecated in a later breaking release
islinear(m::Map) = islinearmap(m)
isaffine(m::Map) = isaffinemap(m)
isconstant(m::Map) = isconstantmap(m)
constant(m::Map) = mapconstant(m)
isidentity(m::Map) = isidentitymap(m)

isreal(::Type{T}) where T = isrealtype(T)

matrix(m::Map) = affinematrix(m)
vector(m::Map) = affinevector(m)

@deprecate convert_domaintype(map::Map, ::Type{T}) where {T} convert_domaintype(T, map)
@deprecate convert_numtype(map::Map{T}, ::Type{U}) where {T,U} convert_numtype(U, map)
@deprecate convert_prectype(map::Map{T}, ::Type{U}) where {T,U} convert_prectype(U, map)
