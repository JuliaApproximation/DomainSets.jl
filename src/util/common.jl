
unitvector(d::Domain{T}, dim) where {N,S,T<:SVector{N,S}} = SVector{N,S}(ntuple(i -> i==dim, N))
function unitvector(d::Domain{T}, dim) where {T<:AbstractVector}
    p = zeros(eltype(T), dimension(d))
    p[dim] = 1
    p
end
unitvector(d::Domain{T}, dim) where {T<:Number} = (@assert dim==1; one(T))

origin(d::Domain{T}) where {T <: StaticTypes} = zero(T)
function origin(d::Domain{T}) where {T <: AbstractVector}
	p = similar(choice(d))
	fill!(p, 0)
	convert(T, p)
end
