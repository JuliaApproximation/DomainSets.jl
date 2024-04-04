"""
    ComplexUnitCircle()
    ComplexUnitCircle{T}()

The unit circle in the complex plane.

See also: [`ComplexUnitDisk`](@ref).
"""
const ComplexUnitCircle{T} = StaticUnitSphere{Complex{T}}

"""
    ComplexUnitDisk()
    ComplexUnitDisk{T}()
    ComplexUnitDisk{T,C}()

The unit disk in the complex plane. The disk is open when `C=:open` and closed
when `C=:closed`.

See also: [`ComplexUnitCircle`](@ref).
"""
const ComplexUnitDisk{T,C} = StaticUnitBall{Complex{T},C}

ComplexUnitCircle() = ComplexUnitCircle{Float64}()
ComplexUnitDisk() = ComplexUnitDisk{Float64}()
ComplexUnitDisk{T}() where {T} = ComplexUnitDisk{T,:closed}()

show(io::IO, d::ComplexUnitCircle{Float64}) = print(io, "ComplexUnitCircle()")
show(io::IO, d::ComplexUnitDisk{Float64,:closed}) = print(io, "ComplexUnitDisk()")
show(io::IO, d::ComplexUnitDisk{Float64,:open}) = print(io, "ComplexUnitDisk()  (open)")
show(io::IO, d::ComplexUnitDisk{T,:closed}) where {T} = print(io, "ComplexUnitDisk{$(T)}()")
show(io::IO, d::ComplexUnitDisk{T,:open}) where {T} = print(io, "ComplexUnitDisk{$(T)}()  (open)")

canonicaldomain(::Parameterization, d::ComplexUnitCircle{T}) where {T} =
    UnitInterval{T}()
mapfrom_canonical(::Parameterization, d::ComplexUnitCircle{T}) where {T} =
    VectorToComplex{T}() ∘ UnitCircleMap{T}()

canonicaldomain(::Isomorphic, d::ComplexUnitCircle{T}) where {T} = UnitCircle{T}()
mapfrom_canonical(::Isomorphic, d::ComplexUnitCircle{T}) where {T} = VectorToComplex{T}()

canonicaldomain(::Parameterization, d::ComplexUnitDisk{T}) where {T} = UnitSquare{T}()
mapfrom_canonical(::Parameterization, d::ComplexUnitDisk{T}) where {T} =
    VectorToComplex{T}() ∘ UnitDiskMap{T}()

canonicaldomain(::Isomorphic, d::ComplexUnitDisk{T}) where {T} = UnitCircle{T}()
mapfrom_canonical(::Isomorphic, d::ComplexUnitDisk{T}) where {T} = VectorToComplex{T}()
