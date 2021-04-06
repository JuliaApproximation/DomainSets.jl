
"An isomorphism is a bijection between types that preserves norms."
abstract type Isomorphism{T,U} <: TypedMap{T,U} end

"Map a length 1 vector `x` to `x[1]`."
struct VectorToNumber{T} <: Isomorphism{SVector{1,T},T}
end
"Map a number `x` to a length 1 vector `[x]`."
struct NumberToVector{T} <: Isomorphism{T,SVector{1,T}}
end

inverse(::VectorToNumber{T}) where {T} = NumberToVector{T}()
inverse(::NumberToVector{T}) where {T} = VectorToNumber{T}()

applymap(::VectorToNumber, x) = x[1]
applymap(::NumberToVector, x) = SVector(x)


"Map a length 2 vector `x` to `x[1] + im*x[2]`."
struct VectorToComplex{T} <: Isomorphism{SVector{2,T},Complex{T}}
end
"Map a complex number `x` to the length 2 vector `[real(x); imag(x)]`."
struct ComplexToVector{T} <: Isomorphism{Complex{T},SVector{2,T}}
end

applymap(::VectorToComplex, x) = x[1] + im*x[2]
applymap(::ComplexToVector, x) = SVector(real(x), imag(x))

inverse(::VectorToComplex{T}) where {T} = ComplexToVector{T}()
inverse(::ComplexToVector{T}) where {T} = VectorToComplex{T}()


"Map a static vector to a tuple."
struct VectorToTuple{N,T} <: Isomorphism{SVector{N,T},NTuple{N,T}}
end
"Map a tuple to a static vector."
struct TupleToVector{N,T} <: Isomorphism{NTuple{N,T},SVector{N,T}}
end

inverse(::VectorToTuple{N,T}) where {N,T} = TupleToVector{N,T}()
inverse(::TupleToVector{N,T}) where {N,T} = VectorToTuple{N,T}()

applymap(::VectorToTuple, x) = tuple(x...)
applymap(::TupleToVector, x) = SVector(x)
