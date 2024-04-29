using DomainSets, CairoMakie, StaticArrays
using DomainSets: Sphere


##
# 2D
##
p = plot((0..1) × (1..2))
plot!(UnitDisk())
plot!(Sphere(2.0, SVector(1.0, 0.5)))
plot!(Sphere(3.0, [1.0, 0.5]))
plot!(Point(SVector(2,3)))
p



##
# 3D
##
p = plot((0..1) × (1..2) × (3..4))
plot!(UnitBall())
plot!(Sphere(2.0, SVector(1.0, 0.5,0.5)))
plot!(Sphere(3.0, [1.0, 0.5,0.5]))
p
