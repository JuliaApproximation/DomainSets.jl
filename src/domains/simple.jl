# simple.jl
# A collection of simple domains.

###########################
# The unit ball and sphere
###########################

"""
The unit ball (of radius 1) in `N` dimensions.
"""
struct UnitBall{N,T} <: EuclideanDomain{N,T}
end

const Disk{T} = UnitBall{2,T}
const Ball{T} = UnitBall{3,T}

indomain(x, ::UnitBall) = norm(x) <= 1

disk(::Type{T} = Float64) where {T} = Disk{T}()
disk(radius::Number) = radius * disk(typeof(radius))
disk(radius::Number, center::AbstractVector) = disk(radius) + center

ball(::Type{T} = Float64) where {T} = Ball{T}()
ball(radius::Number) = radius * ball(typeof(radius))
ball(radius::Number, center::AbstractVector) = ball(radius) + center

show(io::IO, d::UnitBall{N}) where {N} = print(io, "the $(N)-dimensional unit ball")



"""
The unit sphere (of radius 1) in `N` dimensions.
"""
struct UnitSphere{N,T} <: EuclideanDomain{N,T}
end

const Circle{T} = UnitSphere{2,T}
const Sphere{T} = UnitBall{3,T}

indomain(x, ::UnitSphere) = norm(x) == 1

∂(::UnitBall{N,T}) where {N,T} = UnitSphere{N,T}()

circle(::Type{T} = Float64) where {T} = Circle{T}()
circle(radius::Number) = radius * circle(typeof(radius))
circle(radius::Number, center::AbstractVector) = circle(radius) + center

sphere(::Type{T} = Float64) where {T} = Sphere{T}()
sphere(radius::Number) = radius * sphere(typeof(radius))
sphere(radius::Number, center::AbstractVector) = sphere(radius) + center


###########################
# An n-dimensional simplex
###########################

struct UnitSimplex{N,T} <: EuclideanDomain{N,T}
end

indomain(x, ::UnitSimplex) = x .>= 0 && norm(x,1) <= 1

simplex(::Type{Val{N}}, ::Type{T} = Float64) where {T,N} = UnitSimplex{N,T}()


#########################
# An n-dimensional cube
#########################

cube(::Type{Val{N}}, ::Type{T} = Float64) where {N,T} = tensorproduct(UnitInterval{T}(), Val{N})

cube() = cube(Val{3})

rectangle(a, b, c, d) = interval(a,b) ⊗ interval(c,d)

cube(a, b, c, d, e, f) = interval(a,b) ⊗ interval(c,d) ⊗ interval(e,f)


# const Square{T} = UnitCube{2,T}
# const Cube{T} = UnitCube{3,T}}



##############
# A cylinder
##############

cylinder(::Type{T} = Float64) where {T} = disk(T) ⊗ UnitInterval{T}()

cylinder(radius::T, length::T) where {T} = disk(radius) ⊗ Interval(0,length)
