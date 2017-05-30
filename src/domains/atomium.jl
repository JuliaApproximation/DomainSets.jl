# atomium.jl

###
# The atomium: a famous building in Belgium
###

function atomium()
    sphere1 = Ball(0.25)
    spheres = DomainCollection(sphere1)
    push!(spheres, sphere1 + [ 0.6, 0.6, 0.6])
    push!(spheres, sphere1 + [ 0.6, 0.6,-0.6])
    push!(spheres, sphere1 + [ 0.6,-0.6, 0.6])
    push!(spheres, sphere1 + [ 0.6,-0.6,-0.6])
    push!(spheres, sphere1 + [-0.6, 0.6, 0.6])
    push!(spheres, sphere1 + [-0.6, 0.6,-0.6])
    push!(spheres, sphere1 + [-0.6,-0.6, 0.6])
    push!(spheres, sphere1 + [-0.6,-0.6,-0.6])
    cyl1 = Cylinder(0.10, 1.2)
    push!(spheres, cyl1 + [-0.6, 0.6, 0.6]);
    push!(spheres, cyl1 + [-0.6,-0.6, 0.6]);
    push!(spheres, cyl1 + [-0.6, 0.6,-0.6]);
    push!(spheres, cyl1 + [-0.6,-0.6,-0.6]);
    cyl2 = rotate(cyl1, 0.0, 0.0, pi/2.0)
    push!(spheres, cyl2 + [ 0.6, -0.6, 0.6])
    push!(spheres, cyl2 + [-0.6, -0.6, 0.6])
    push!(spheres, cyl2 + [ 0.6, -0.6,-0.6])
    push!(spheres, cyl2 + [-0.6, -0.6,-0.6])
    cyl2b = rotate(cyl1, 0.0, pi/2.0, 0.0)
    push!(spheres, cyl2b + [ 0.6,  0.6, 0.6])
    push!(spheres, cyl2b + [-0.6,  0.6, 0.6])
    push!(spheres, cyl2b + [ 0.6, -0.6, 0.6])
    push!(spheres, cyl2b + [-0.6, -0.6, 0.6])
    cyl3 = Cylinder(0.10, 1.2*sqrt(3))
    cyl3 = rotate(cyl3, 0.0, asin(1/sqrt(3)), 0.0)
    cyl3 = rotate(cyl3, 0.0, 0.0, pi/4)
    push!(spheres, cyl3 + [ -0.6, -0.6, +0.6])
    cyl4 = Cylinder(0.10, 1.2*sqrt(3))
    cyl4 = rotate(cyl4, 0.0, -asin(1/sqrt(3)), 0.0)
    cyl4 = rotate(cyl4, 0.0, 0.0, pi/4)
    push!(spheres, cyl4 + [ -0.6, -0.6, -0.6])
    cyl5 = Cylinder(0.10, 1.2*sqrt(3))
    cyl5 = rotate(cyl5, 0.0, asin(1/sqrt(3)), 0.0)
    cyl5 = rotate(cyl5, 0.0, 0.0, -pi/4)
    push!(spheres, cyl5 + [ -0.6, +0.6, +0.6])
    cyl6 = Cylinder(0.10, 1.2*sqrt(3))
    cyl6 = rotate(cyl6, 0.0, -asin(1/sqrt(3)), 0.0)
    cyl6 = rotate(cyl6, 0.0, 0.0, -pi/4)
    push!(spheres, cyl6 + [ -0.6, +0.6, -0.6])
    spheres.box = unitbox3
    atomium = spheres
end
