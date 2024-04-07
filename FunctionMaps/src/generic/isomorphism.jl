
"""
    Isomorphism{T,U} <: TypedMap{T,U}

An isomorphism is a bijection between types that preserves norms.
"""
abstract type Isomorphism{T,U} <: TypedMap{T,U} end

show(io::IO, m::Isomorphism{T,U}) where {T,U} = print(io, "x : $(T) -> x : $(U)")
map_object_parentheses(m::Isomorphism) = true

struct VectorToNumber{T} <: Isomorphism{SVector{1,T},T}
end
struct NumberToVector{T} <: Isomorphism{T,SVector{1,T}}
end

"""
    VectorToNumber()
    VectorToNumber{T}()

Map a length 1 vector `x` to the number `x[1]`.

See also: [`NumberToVector`](@ref).
"""
VectorToNumber() = VectorToNumber{Float64}()

"""
    NumberToVector()
    NumberToVector{T}()

Map a number `x` to the length 1 vector `[x]`.

See also: [`VectorToNumber`](@ref).
"""
NumberToVector() = NumberToVector{Float64}()

mapsize(::VectorToNumber) = (1,1)
mapsize(::NumberToVector) = (1,)

inverse(::VectorToNumber{T}) where {T} = NumberToVector{T}()
inverse(::NumberToVector{T}) where {T} = VectorToNumber{T}()
inverse(m::VectorToNumber, x) = inverse(m)(x)
inverse(m::NumberToVector, x) = inverse(m)(x)

applymap(::VectorToNumber, x) = x[1]
applymap(::NumberToVector, x) = SVector(x)

jacobian(::VectorToNumber{T}, x) where {T} = transpose(SVector(one(T)))
jacobian(::VectorToNumber{T}) where {T} = ConstantMap{SVector{1,T}}(transpose(SVector(one(T))))

jacobian(::NumberToVector{T}, x) where {T} = SVector(one(T))
jacobian(::NumberToVector{T}) where {T} = ConstantMap{T}(SVector(one(T)))


struct VectorToComplex{T} <: Isomorphism{SVector{2,T},Complex{T}}
end
struct ComplexToVector{T} <: Isomorphism{Complex{T},SVector{2,T}}
end

"""
    VectorToComplex()
    VectorToComplex{T}()

Map a length 2 vector ``[a;b]`` to the complex number ``a+bi``.

See also: [`ComplexToVector`](@ref).
"""
VectorToComplex() = VectorToComplex{Float64}()

"""
    ComplexToVector()
    ComplexToVector{T}()

Map a complex number ``a+bi`` to the length 2 vector ``[a; b]``.

See also: [`VectorToComplex`](@ref).
"""
ComplexToVector() = ComplexToVector{Float64}()

mapsize(::VectorToComplex) = (1,2)
mapsize(::ComplexToVector) = (2,)

applymap(::VectorToComplex, x) = x[1] + im*x[2]
applymap(::ComplexToVector, x) = SVector(real(x), imag(x))

inverse(::VectorToComplex{T}) where {T} = ComplexToVector{T}()
inverse(::ComplexToVector{T}) where {T} = VectorToComplex{T}()
inverse(m::VectorToComplex, x) = inverse(m)(x)
inverse(m::ComplexToVector, x) = inverse(m)(x)

jacobian(::VectorToComplex{T}, x) where {T} = transpose(SVector(one(T),one(T)*im))
jacobian(::VectorToComplex{T}) where {T} = ConstantMap{SVector{2,T}}(transpose(SVector(one(T),one(T)*im)))


"Map a static vector to a tuple."
struct VectorToTuple{N,T} <: Isomorphism{SVector{N,T},NTuple{N,T}}
end
"Map a tuple to a static vector."
struct TupleToVector{N,T} <: Isomorphism{NTuple{N,T},SVector{N,T}}
end

inverse(::VectorToTuple{N,T}) where {N,T} = TupleToVector{N,T}()
inverse(::TupleToVector{N,T}) where {N,T} = VectorToTuple{N,T}()
inverse(m::VectorToTuple, x) = inverse(m)(x)
inverse(m::TupleToVector, x) = inverse(m)(x)

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
inverse(m::NestedToFlat, x) = inverse(m)(x)
inverse(m::FlatToNested, x) = inverse(m)(x)

applymap(::NestedToFlat{N,T,U,DIM}, x) where {N,T,U,DIM} = convert_tocartesian(x, Val{DIM}())
applymap(::FlatToNested{N,T,U,DIM}, x) where {N,T,U,DIM} = convert_fromcartesian(x, Val{DIM}())
