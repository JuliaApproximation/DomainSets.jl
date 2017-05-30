# simple.jl

################################################################################
### The unit ball in N dimensions
################################################################################

struct UnitBall{N} <: Domain{N}
end

indomain(x, ::UnitBall{1}) = -1 <= x <= 1
indomain(x, ::UnitBall{2}) = x[1]^2+x[2]^2 <= 1
indomain(x, ::UnitBall{3}) = x[1]^2+x[2]^2+x[3]^2 <= 1

indomain{N}(x, ::UnitBall{N}) = sum(map(t->t^2, x)) <= 1

boundingbox{N}(::UnitBall{N}) = BBox{N,Int}(-ones(SVector{N,Int}), ones(SVector{N,Int}))

Disk() = UnitBall{2}()
Disk(radius) = radius * Disk()
Disk(radius, center) = radius * Disk() + center

show(io::IO, d::UnitBall) = print(io, "the $(ndims(d))-dimensional unit ball")

const unitdisk = Disk()



################################################################################
### A 3D ball
################################################################################

struct Ball{S,T} <: Domain{3}
    radius    ::  S
    center    ::  SVector{3,T}

    Ball{S,T}(radius = one(S), center = zeros(SVector{3,T})) where {S,T} = new(radius, center)
end

Ball() = Ball{Int,Float64}()
Ball{T}(::Type{T}) = Ball{T,T}()

Ball{T}(radius::T) = Ball{T,T}(radius)
Ball{S,T}(radius::S, center::SVector{3,T}) = Ball{S,T}(radius, center)
Ball(radius, center::AbstractVector) = Ball(radius, SVector{3}(center))


indomain(x, s::Ball) = (x[1]-s.center[1])^2 + (x[2]-s.center[2])^2 + (x[3]-s.center[3])^2 <= s.radius^2

## Arithmetic operations

(+)(s::Ball, x::SVector{3}) = Ball(s.radius, s.center+x)

(*)(s::Ball, x::Number) = Ball(s.radius * x, s.center * x)


boundingbox(c::Ball) = BBox((c.center[1]-c.radius,c.center[2]-c.radius,c.center[3]-c.radius),(c.center[1]+c.radius,c.center[2]+c.radius,c.center[3]+c.radius))

show(io::IO, s::Ball) = print(io, "a ball of radius ", s.radius, " centered at ", s.center)

const unitball = Ball()





################################################################################
### An n-dimensional cube
################################################################################

Cube() = Cube(Val{1})

Cube(::Type{Val{1}}) = Interval()
Cube{N}(::Type{Val{N}}) = Interval() ⊗ Cube(Val{N-1})
Cube(n::Int) = Cube(Val{n})

Cube(left::Number, right::Number) = Interval(left, right)
Cube(left, right) = tensorproduct(map(Interval, left, right)...)

rectangle(a, b, c, d) = Interval(a,b) ⊗ Interval(c,d)

cube(a, b, c, d, e, f) = Interval(a,b) ⊗ Interval(c,d) ⊗ Interval(e,f)


const unitsquare = Cube(Val{2})
const unitcube = Cube(Val{3})



################################################################################
### A cylinder
################################################################################


cylinder(radius = 1, length = 1) = Disk(radius) ⊗ Interval(0,length)



##########################################################################
### Assorted Domains
##########################################################################


function randomcircles(n, radius = 0.3)
    list = Domain[Disk(radius, (2*rand(2)-1)*0.8) for i=1:n]
    DomainCollection{2}(list)
end
