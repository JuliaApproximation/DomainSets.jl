
"A lazy Jacobian `J` stores a map `m` and returns `J(x) = jacobian(m, x)`."
struct LazyJacobian{T,M} <: SimpleLazyMap{T}
	map	::	M
end

LazyJacobian(m::Map{T}) where {T} = LazyJacobian{T,typeof(m)}(m)
applymap(m::LazyJacobian, x) = jacobian(supermap(m), x)

size(m::LazyJacobian) = size(supermap(m))

Display.displaystencil(m::LazyJacobian) = ["LazyJacobian(", supermap(m), ")"]
show(io::IO, mime::MIME"text/plain", m::LazyJacobian) = composite_show(io, mime, m)


"""
    jacobian(m::AbstractMap[, x])

Return the jacobian map. The two-argument version evaluates the jacobian
at a point `x`.
"""
jacobian(m::AbstractMap) = LazyJacobian(m)
jacobian!(y, m::AbstractMap, x) = y .= jacobian(m, x)


"A `DeterminantMap` returns the determinant of the result of a given map."
struct DeterminantMap{T,M} <: SimpleLazyMap{T}
    map ::  M
end

DeterminantMap(m::Map{T}) where {T} = DeterminantMap{T}(m)
DeterminantMap{T}(m) where {T} = DeterminantMap{T,typeof(m)}(m)
DeterminantMap{T}(m::Map{T}) where {T} = DeterminantMap{T,typeof(m)}(m)
DeterminantMap{T}(m::Map{S}) where {S,T} = DeterminantMap{T}(convert(Map{T}, m))

determinantmap(m::AbstractMap) = DeterminantMap(m)

applymap(m::DeterminantMap, x) = det(supermap(m)(x))


"""
    jacdet(m::AbstractMap[, x])

Return the determinant of the jacobian as a map. The two-argument version
evaluates the jacobian determinant at a point `x`.
"""
jacdet(m::AbstractMap) = determinantmap(jacobian(m))
jacdet(m::AbstractMap, x) = det(jacobian(m, x))


"An `AbsMap` returns the absolute value of the result of a given map."
struct AbsMap{T,M} <: SimpleLazyMap{T}
    map ::  M
end

AbsMap(m::Map{T}) where {T} = AbsMap{T}(m)
AbsMap{T}(m) where {T} = AbsMap{T,typeof(m)}(m)
AbsMap{T}(m::Map{T}) where {T} = AbsMap{T,typeof(m)}(m)
AbsMap{T}(m::Map{S}) where {S,T} = AbsMap{T}(convert(Map{T}, m))

absmap(m::AbstractMap) = AbsMap(m)

applymap(m::AbsMap, x) = abs(supermap(m)(x))


"A lazy volume element evaluates to `diffvolume(m, x)` on the fly."
struct LazyDiffVolume{T,M} <: SimpleLazyMap{T}
	map	::	M
end

LazyDiffVolume(m::Map{T}) where {T} = LazyDiffVolume{T,typeof(m)}(m)
applymap(m::LazyDiffVolume, x) = diffvolume(supermap(m), x)

Display.displaystencil(m::LazyDiffVolume) = ["LazyDiffVolume(", supermap(m), ")"]
show(io::IO, mime::MIME"text/plain", m::LazyDiffVolume) = composite_show(io, mime, m)

"""
    diffvolume(m::AbstractMap[, x])

Compute the differential volume (at a point `x`). If `J` is the Jacobian matrix,
possibly rectangular, then the differential volume is `sqrt(det(J'*J))`.

If the map is square, then the differential volume is the absolute value of the
Jacobian determinant.
"""
function diffvolume(m::AbstractMap)
	if issquare(m)
		if size(m,1) > 1
			absmap(jacdet(m))
		else
			jacdet(m)
		end
	else
		LazyDiffVolume(m)
	end
end

function diffvolume(m::AbstractMap, x)
	if issquare(m)
		if size(m,1) > 1
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
to_matrix(::Type{SVector{N,T}}, A::Number) where {N,T} = A * one(SMatrix{N,N,T})
to_matrix(::Type{SVector{N,T}}, A::UniformScaling) where {N,T} = A.λ * one(SMatrix{N,N,T})
to_matrix(::Type{T}, A::Number) where {T<:AbstractVector} = A * I
to_matrix(::Type{T}, A::UniformScaling) where {T<:Number} = one(T)
to_matrix(::Type{T}, A::UniformScaling) where {T<:AbstractVector} = A

to_matrix(::Type{T}, A, b) where {T} = A
to_matrix(::Type{T}, A::AbstractMatrix, b) where {T} = A
to_matrix(::Type{T}, A::Number, b::Number) where {T<:Number} = A
to_matrix(::Type{T}, A::UniformScaling, b::Number) where {T<:Number} = A.λ
to_matrix(::Type{SVector{N,T}}, A::NumberLike, b::SVector{N,T}) where {N,T} = A * one(SMatrix{N,N,T})
to_matrix(::Type{T}, A::NumberLike, b::AbstractVector) where {S,T<:AbstractVector{S}} =
    A * Array{S,2}(I, length(b), length(b))

"""
    to_vector(::Type{T}, A[, b])

Convert the `b` in the affine map `A*x` or `A*x+b` with domaintype `T` to a vector.
"""
to_vector(::Type{T}, A) where {T} = zero(T)
to_vector(::Type{T}, A::SVector{M,S}) where {T,M,S} = zero(SVector{M,S})
to_vector(::Type{T}, A::SMatrix{M,N,S}) where {T<:AbstractVector,M,N,S} = zero(SVector{M,S})
to_vector(::Type{T}, A::AbstractArray) where {T<:AbstractVector} = zeros(eltype(T),size(A,1))
to_vector(::Type{T}, A, b) where {T} = b


"Return a zero matrix of the same size as the map."
zeromatrix(m::Map) = zeromatrix(m, domaintype(m), codomaintype(m))
zeromatrix(m, ::Type{T}, ::Type{U}) where {T,U} = zeros(numtype(T),size(m))
zeromatrix(m, ::Type{T}, ::Type{U}) where {T<:Number,U<:Number} = zero(promote_type(T,U))
zeromatrix(m, ::Type{T}, ::Type{U}) where {N,T<:StaticVector{N},U<:Number} =
	transpose(zero(SVector{N,eltype(T)}))
zeromatrix(m, ::Type{T}, ::Type{U}) where {N,M,T<:StaticVector{N},U<:StaticVector{M}} =
	zero(SMatrix{M,N,promote_type(eltype(T),eltype(U))})
zeromatrix(m, ::Type{T}, ::Type{U}) where {T<:Number,M,U<:StaticVector{M}} =
	zero(SVector{M,promote_type(T,eltype(U))})
zeromatrix(m, ::Type{T}, ::Type{U}) where {T<:AbstractVector,U<:Number} =
	transpose(zeros(promote_type(eltype(T),U), size(m,1)))


"Return a zero vector of the same size as the codomain of the map."
zerovector(m::Map) = zerovector(m, codomaintype(m))
zerovector(m::Map, ::Type{U}) where {U} = zero(U)
zerovector(m::Map, ::Type{StaticVector{M,T}}) where {M,T} = zero(SVector{M,T})
# If the output type is a vector, the map itself should store the size information.
zerovector(m::Map, ::Type{<:AbstractVector{T}}) where {T} = zeros(T, size(m,1))

"Return an identity matrix with the same size as the map."
identitymatrix(m::Map) = identitymatrix(m, codomaintype(m))
identitymatrix(m::Map, ::Type{T}) where {T} = one(T)
identitymatrix(m::Map, ::Type{<:StaticVector{N,T}}) where {N,T} = one(SMatrix{N,N,T})
identitymatrix(m::Map, ::Type{<:AbstractVector{T}}) where {T} = Diagonal{T}(ones(size(m,1)))
