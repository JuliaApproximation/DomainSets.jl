
# This file contains a few specific maps, included mainly to test
# functionality on non-trivial maps

"""
A Cartesion to Polar map. First dimension is interpreted as radial distance,
second as an angle.
The unit circle is mapped to the square `[-1,1]x[-1,1]`.
"""
struct CartToPolarMap{T} <: Map{SVector{2,T}}
end

CartToPolarMap() = CartToPolarMap{Float64}()

mapsize(m::CartToPolarMap) = (2,2)

applymap(map::CartToPolarMap{T}, x) where {T} =
    SVector{2,T}(sqrt(x[1]^2+x[2]^2)*2-1, atan(x[2],x[1])/pi)

function jacobian(m::CartToPolarMap{T}, x) where {T}
    d = sqrt(x[1]^2+x[2]^2)
    r = 1/(1+(x[2]/x[1])^2)
    SMatrix{2,2,T}(2*x[1]/d, 1/pi*r*(-x[2]/x[1]^2), 2*x[2]/d, 1/pi*r*1/x[1])
end

inverse(m::CartToPolarMap{T}) where {T} = PolarToCartMap{T}()
inverse(m::CartToPolarMap, x) = inverse(m)(x)

isreal(m::CartToPolarMap) = true

convert(::Type{Map{SVector{2,T}}}, ::CartToPolarMap) where {T} = CartToPolarMap{T}()

isequalmap(m1::CartToPolarMap, m2::CartToPolarMap) = true


"""
A Polar to Cartesian map. The angle is mapped to the second dimension,
radius to the first.
The square `[-1,1]x[-1,1]` is mapped to the unit circle.
"""
struct PolarToCartMap{T} <: Map{SVector{2,T}}
end

PolarToCartMap() = PolarToCartMap{Float64}()

mapsize(m::PolarToCartMap) = (2,2)

applymap(map::PolarToCartMap{T}, x) where {T} = SVector{2,T}((x[1]+1)/2*cos(pi*x[2]), (x[1]+1)/2*sin(pi*x[2]))

jacobian(m::PolarToCartMap{T}, x) where {T} =
    SMatrix{2,2,T}(cos(pi*x[2])/2, sin(pi*x[2])/2, -pi*(x[1]+1)/2*sin(pi*x[2]), pi*(x[1]+1)/2*cos(pi*x[2]))

inverse(m::PolarToCartMap{T}) where {T} = CartToPolarMap{T}()
inverse(m::PolarToCartMap, x) = inverse(m)(x)

isreal(m::PolarToCartMap) = true

convert(::Type{Map{SVector{2,T}}}, ::PolarToCartMap) where {T} = PolarToCartMap{T}()

isequalmap(m1::PolarToCartMap, m2::PolarToCartMap) = true


"""
The map `[cos(2πt), sin(2πt)]` from `[0,1]` to the unit circle in `ℝ^2`.
"""
struct UnitCircleMap{T} <: Map{T} end

UnitCircleMap() = UnitCircleMap{Float64}()

mapsize(m::UnitCircleMap) = (2,)

applymap(m::UnitCircleMap{T}, t) where {T} = SVector(cos(2*T(pi)*t), sin(2*T(pi)*t))

function jacobian(m::UnitCircleMap{T}, t) where {T}
    a = 2*T(pi)
    SVector(-a*sin(a*t), a*cos(a*t))
end

# we know that the differential volume is 2*pi
diffvolume(m::UnitCircleMap{T}) where {T} = ConstantMap{T}(2*T(pi))

"""
`AngleMap` is a left inverse of `UnitCircleMap`. A 2D vector `x` is projected onto
the intersection point with the unit circle of the line connecting `x` to the
origin. The angle of this point, scaled to the interval `[0,1)`, is the result.
"""
struct AngleMap{T} <: Map{SVector{2,T}}
end

AngleMap() = AngleMap{Float64}()

function applymap(m::AngleMap{T}, x) where {T}
    twopi = 2*convert(T, pi)
    θ = atan(x[2],x[1])
    if θ < 0
        # atan2 returns an angle in (-π,π], convert to [0,2π) using periodicity.
        θ += twopi
    end
    # And divide by 2π to scale to [0,1)
    θ / twopi
end

mapsize(m::AngleMap) = (1,2)

function jacobian(m::AngleMap{T}, t) where {T}
    x = t[1]; y = t[2]
    twopi = 2*convert(T, pi)
    v = SVector(one(T) / (1+(y/x)^2) * (-y/x^2) * 1/twopi,
            one(T) / (1+(y/x)^2) * one(T)/x * 1/twopi)
    transpose(v)
end


leftinverse(m::UnitCircleMap{T}) where {T} = AngleMap{T}()
leftinverse(m::UnitCircleMap, x) = leftinverse(m)(x)
rightinverse(m::AngleMap{T}) where {T} = UnitCircleMap{T}()
rightinverse(m::AngleMap, x) = rightinverse(m)(x)


"""
The map `r*[cos(2πt), sin(2πt)]` from `[0,1]^2` to the unit disk in `ℝ^2`.
"""
struct UnitDiskMap{T} <: Map{SVector{2,T}} end

UnitDiskMap() = UnitDiskMap{Float64}()

mapsize(m::UnitDiskMap) = (2,2)

applymap(m::UnitDiskMap{T}, x) where {T} =
    SVector(x[1]*cos(2*T(pi)*x[2]), x[1]*sin(2*T(pi)*x[2]))

function jacobian(m::UnitDiskMap{T}, x) where {T}
    a = 2*T(pi)
    SMatrix{2,2}(cos(a*x[2]), sin(a*x[2]), -a*x[1]*sin(a*x[2]), a*x[1]*cos(a*x[2]))
end

# we know that the differential volume is 2*r*pi
diffvolume(m::UnitDiskMap, x) = 2*x[1]*pi
