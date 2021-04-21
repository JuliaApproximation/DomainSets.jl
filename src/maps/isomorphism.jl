
"An isomorphism is a bijection between types that preserves norms."
abstract type Isomorphism{T,U} <: TypedMap{T,U} end

show(io::IO, m::Isomorphism{T,U}) where {T,U} = print(io, "x : $(T) -> x : $(U)")
Display.object_parentheses(m::Isomorphism) = true

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


"Map a nested vector or tuple to a flat vector."
struct NestedToFlat{N,T,U,DIM} <: Isomorphism{U,SVector{N,T}}
end
"Map a flattened vector to a nested one."
struct FlatToNested{N,T,U,DIM} <: Isomorphism{SVector{N,T},U}
end

inverse(::NestedToFlat{N,T,U,DIM}) where {N,T,U,DIM} = FlatToNested{N,T,U,DIM}()
inverse(::FlatToNested{N,T,U,DIM}) where {N,T,U,DIM} = NestedToFlat{N,T,U,DIM}()

applymap(::NestedToFlat{N,T,U,DIM}, x) where {N,T,U,DIM} = convert_tocartesian(x, Val{DIM}())
applymap(::FlatToNested{N,T,U,DIM}, x) where {N,T,U,DIM} = convert_fromcartesian(x, Val{DIM}())
