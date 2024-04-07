
"A lazy Jacobian `J` stores a map `m` and returns `J(x) = jacobian(m, x)`."
struct LazyJacobian{T,M} <: SimpleLazyMap{T}
	map	::	M
end

LazyJacobian(m) = LazyJacobian{domaintype(m),typeof(m)}(m)
applymap(m::LazyJacobian, x) = jacobian(supermap(m), x)

mapsize(m::LazyJacobian) = mapsize(supermap(m))

Display.displaystencil(m::LazyJacobian) = ["LazyJacobian(", supermap(m), ")"]
show(io::IO, mime::MIME"text/plain", m::LazyJacobian) = composite_show(io, mime, m)


"""
    jacobian(m[, x])

Return the jacobian map. The two-argument version evaluates the jacobian
at a point `x`.
"""
jacobian(m) = LazyJacobian(m)
jacobian!(y, m, x) = y .= jacobian(m, x)


"A `DeterminantMap` returns the determinant of the result of a given map."
struct DeterminantMap{T,M} <: SimpleLazyMap{T}
    map ::  M
end

DeterminantMap(m) = DeterminantMap{domaintype(m)}(m)
DeterminantMap{T}(m) where {T} = DeterminantMap{T,typeof(m)}(m)
DeterminantMap{T}(m::Map{T}) where {T} = DeterminantMap{T,typeof(m)}(m)
DeterminantMap{T}(m::Map{S}) where {S,T} = DeterminantMap{T}(convert(Map{T}, m))

determinantmap(m) = DeterminantMap(m)

applymap(m::DeterminantMap, x) = det(supermap(m)(x))


"""
    jacdet(m[, x])

Return the determinant of the jacobian as a map. The two-argument version
evaluates the jacobian determinant at a point `x`.
"""
jacdet(m) = determinantmap(jacobian(m))
jacdet(m, x) = det(jacobian(m, x))


"An `AbsMap` returns the absolute value of the result of a given map."
struct AbsMap{T,M} <: SimpleLazyMap{T}
    map ::  M
end

AbsMap(m) = AbsMap{domaintype(m)}(m)
AbsMap{T}(m) where {T} = AbsMap{T,typeof(m)}(m)
AbsMap{T}(m::Map{T}) where {T} = AbsMap{T,typeof(m)}(m)
AbsMap{T}(m::Map{S}) where {S,T} = AbsMap{T}(convert(Map{T}, m))

absmap(m) = AbsMap(m)

applymap(m::AbsMap, x) = abs(supermap(m)(x))


"A lazy volume element evaluates to `diffvolume(m, x)` on the fly."
struct LazyDiffVolume{T,M} <: SimpleLazyMap{T}
	map	::	M
end

LazyDiffVolume(m) = LazyDiffVolume{domaintype(m)}(m)
LazyDiffVolume{T}(m) where {T} = LazyDiffVolume{T,typeof(m)}(m)
LazyDiffVolume{T}(m::Map{T}) where {T} = LazyDiffVolume{T,typeof(m)}(m)
LazyDiffVolume{T}(m::Map{S}) where {S,T} = LazyDiffVolume{T}(convert(Map{T}, m))

applymap(m::LazyDiffVolume, x) = diffvolume(supermap(m), x)

Display.displaystencil(m::LazyDiffVolume) = ["LazyDiffVolume(", supermap(m), ")"]
show(io::IO, mime::MIME"text/plain", m::LazyDiffVolume) = composite_show(io, mime, m)

"""
    diffvolume(m[, x])

Compute the differential volume (at a point `x`). If `J` is the Jacobian matrix,
possibly rectangular, then the differential volume is `sqrt(det(J'*J))`.

If the map is square, then the differential volume is the absolute value of the
Jacobian determinant.
"""
function diffvolume(m)
	if issquaremap(m)
		if mapsize(m,1) > 1
			absmap(jacdet(m))
		else
			jacdet(m)
		end
	else
		LazyDiffVolume(m)
	end
end

function diffvolume(m, x)
	if issquaremap(m)
		if mapsize(m,1) > 1
			abs(jacdet(m, x))
		else
			jacdet(m, x)
		end
	else
		jac = jacobian(m, x)
		sqrt(det(adjoint(jac)*jac))
	end
end



const NumberLike = Union{Number,UniformScaling}

"""
    to_matrix(::Type{T}, A[, b])

Convert the `A` in the affine map `A*x` or `A*x+b` with domaintype `T` to a matrix.
"""
to_matrix(::Type{T}, A) where {T} = A
to_matrix(::Type{T}, A::AbstractMatrix) where {T} = A
to_matrix(::Type{T}, A::NumberLike) where {T<:Number} = A
to_matrix(::Type{<:StaticVector{N,T}}, A::Number) where {N,T} = A * one(SMatrix{N,N,T})
to_matrix(::Type{<:StaticVector{N,T}}, A::UniformScaling) where {N,T} = A.λ * one(SMatrix{N,N,T})
to_matrix(::Type{T}, A::Number) where {T<:AbstractVector} = A * I
to_matrix(::Type{T}, A::UniformScaling) where {T<:Number} = one(T)
to_matrix(::Type{T}, A::UniformScaling) where {T<:AbstractVector} = A

to_matrix(::Type{T}, A, b) where {T} = A
to_matrix(::Type{T}, A::AbstractMatrix, b) where {T} = A
to_matrix(::Type{T}, A::Number, b::Number) where {T<:Number} = A
to_matrix(::Type{T}, A::UniformScaling{S}, b::Number) where {S,T<:Number} =
	convert(promote_type(S,T,typeof(b)), A.λ)
to_matrix(::Type{<:StaticVector{N,T}}, A::NumberLike, b::StaticVector{N,T}) where {N,T} = A * one(SMatrix{N,N,T})
to_matrix(::Type{T}, A::NumberLike, b::AbstractVector) where {S,T<:AbstractVector{S}} =
    A * Array{S,2}(I, length(b), length(b))

"""
    to_vector(::Type{T}, A[, b])

Convert the `b` in the affine map `A*x` or `A*x+b` with domaintype `T` to a vector.
"""
to_vector(::Type{T}, A) where {T} = zero(T)
to_vector(::Type{T}, A::StaticVector{M,S}) where {T,M,S} = zero(SVector{M,S})
to_vector(::Type{T}, A::StaticMatrix{M,N,S}) where {T,M,N,S} = zero(SVector{M,S})
to_vector(::Type{T}, A::AbstractArray) where {T} = zeros(eltype(T),size(A,1))
to_vector(::Type{T}, A, b) where {T} = b


"Return a zero matrix of the same size as the map."
zeromatrix(m) = zeromatrix(m, domaintype(m), codomaintype(m))
zeromatrix(m, ::Type{T}, ::Type{U}) where {T,U} = zeros(numtype(T),mapsize(m))
zeromatrix(m, ::Type{T}, ::Type{U}) where {T<:Number,U<:Number} = zero(promote_type(T,U))
zeromatrix(m, ::Type{T}, ::Type{U}) where {N,T<:StaticVector{N},U<:Number} =
	transpose(zero(SVector{N,eltype(T)}))
zeromatrix(m, ::Type{T}, ::Type{U}) where {N,M,T<:StaticVector{N},U<:StaticVector{M}} =
	zero(SMatrix{M,N,promote_type(eltype(T),eltype(U))})
zeromatrix(m, ::Type{T}, ::Type{U}) where {T<:Number,M,U<:StaticVector{M}} =
	zero(SVector{M,promote_type(T,eltype(U))})
zeromatrix(m, ::Type{T}, ::Type{U}) where {T<:AbstractVector,U<:Number} =
	transpose(zeros(promote_type(eltype(T),U), mapsize(m,1)))


"Return a zero vector of the same size as the codomain of the map."
zerovector(m) = zerovector(m, codomaintype(m))
zerovector(m, ::Type{U}) where {U} = zero(U)
zerovector(m, ::Type{<:StaticVector{M,T}}) where {M,T} = zero(SVector{M,T})
# If the output type is a vector, the map itself should store the size information.
zerovector(m, ::Type{<:AbstractVector{T}}) where {T} = zeros(T, mapsize(m,1))

"Return an identity matrix with the same size as the map."
identitymatrix(m) = identitymatrix(m, codomaintype(m))
identitymatrix(m, ::Type{T}) where {T} = one(T)
identitymatrix(m, ::Type{<:StaticVector{N,T}}) where {N,T} = one(SMatrix{N,N,T})
identitymatrix(m, ::Type{<:AbstractVector{T}}) where {T} = Diagonal{T}(ones(mapsize(m,1)))
