# basic_spaces.jl

"""
A univariate space is the simplest possible space. It only has an element type
`T`, and all possible values of `T` belong to the space.
"""
struct UnivariateSpace{T} <: GeometricSpace{T}
end

UnivariateSpace() = UnivariateSpace{Float64}()


widen(space::UnivariateSpace) = UnivariateSpace{widen(eltype(space))}()

# The symbol is \BbbZ
const IntegerSpace{T <: Integer} = UnivariateSpace{T}
const ℤ = IntegerSpace{Int}()

# The symbol is \BbbR
const RealSpace{T <: AbstractFloat} = UnivariateSpace{T}
const ℝ = RealSpace{Float64}()

# The symbol is \BbbC
const ComplexPlane{T <: AbstractFloat} = UnivariateSpace{Complex{T}}
const ℂ = ComplexPlane{Float64}()

# All univariate spaces are embedded to the same space with a broader T
embedded{S <: UnivariateSpace}(s1::S, s2::S) = promote_type(eltype(s1),eltype(s2)) == eltype(s2)

embedded{T1,T2}(s1::IntegerSpace{T1}, s2::RealSpace{T2}) = promote_type(T1,T2) == T2

embedded{T1,T2}(s1::RealSpace{T1}, s2::ComplexPlane{T2}) = promote_type(T1,T2) == T2
