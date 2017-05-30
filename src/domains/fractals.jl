# fractals.jl

###############################################################################
## The Mandelbrot set
###############################################################################

struct Mandelbrot{T} <: Domain{2}
    maxiter     ::  Int
    threshold   ::  T
    maskcache   ::  Dict
    box         ::  BBox2{T}
end

function Mandelbrot(maxiter = 1000, threshold = 1000.0)
    box = BBox(-1.0, 0.35, -0.65, 0.65)
    M1 = 136
    M2 = 200
    mask1 = computemandelbrotgrid(equispaced_grid(box, (M1,M1)), maxiter, threshold)
    mask2 = computemandelbrotgrid(equispaced_grid(box, (M2,M2)), maxiter, threshold)
    cache = Dict{Int,Array{Bool,2}}()
    cache[M1] = mask1
    cache[M2] = mask2
    Mandelbrot(maxiter, threshold, cache, box)
end


function mandelbrotiteration(x, maxiter, threshold)
    c = 2*(x[1]+1im*x[2])
    T = typeof(c)
    z = zero(T)
    iter = 0
    while abs(z) < threshold && iter < maxiter
        z = z^2 + c
        iter += 1
    end
    abs(z) < threshold
end

function computemandelbrotgrid(grid, maxiter, threshold)
    mask = zeros(Bool, size(grid))
    for i_2 = 1:size(grid, 2)
        for i_1 = 1:size(grid, 1)
            mask[i_1,i_2] = mandelbrotiteration(grid[i_1,i_2], maxiter, threshold)
        end
    end
    mask
end

function indomain_grid(grid, m::Mandelbrot)
    if (left(grid) ≈ left(m.box)) && (right(grid) ≈ right(m.box))
        if haskey(m.maskcache, size(grid,1))
            mask = m.maskcache[size(grid,1)]
        else # compute mask and cache it
            mask = computemandelbrotgrid(grid, m.maxiter, m.threshold)
            m.maskcache[size(grid,1)] = mask
        end
    else # Don't cache if the grid doesn't match the bounding box
        mask = computemandelbrotgrid(grid, m.maxiter, m.threshold)
    end
    mask
end

function isapprox{T}(t::Tuple{T,T}, v::SVector{2,T})
    return t[1]≈v[1] && t[2]≈v[2]
end
function isapprox{T}(v::SVector{2,T}, t::Tuple{T,T})
    return t[1]≈v[1] && t[2]≈v[2]
end
boundingbox(m::Mandelbrot) = m.box

indomain(x::SVector{2}, m::Mandelbrot) = mandelbrotiteration(x, m.maxiter, m.threshold)

show(io::IO, m::Mandelbrot) = print(io, "The Mandelbrot set")


################
## Julia sets
################

struct JuliaSet{T} <: Domain{2}
    c           ::  Complex{T}
    maxiter     ::  Int
    maskcache   ::  Dict
    box         ::  BBox2{T}
end

function JuliaSet(c = -0.122565+0.744866im, maxiter = 1000)
    box = BBox(-0.2, 1.2, -0.4, 0.4)

    mask1 = computejuliasetgrid(equispaced_grid(box, (100,100)), c, maxiter)
    mask2 = computejuliasetgrid(equispaced_grid(box, (200,200)), c, maxiter)
    cache = Dict{Int,Array{Bool,2}}()
    cache[100] = mask1
    cache[200] = mask2
    JuliaSet(c, maxiter, cache, box)
end


function juliasetiteration(x, c, maxiter)
    gamma = 1 + sqrt(1-4*c)
    z = x[1] + 1im*x[2]
    for i = 1:maxiter
        z = gamma*z*(1-z)
    end
    abs(z) < 1000
end

function computejuliasetgrid(grid, c, maxiter)
    m = size(grid)
    mask = zeros(Bool, m)
    for i_2 = 1:m[2]
        for i_1 = 1:m[1]
            a,b = grid[i_1,i_2]
            mask[i_1,i_2] = juliasetiteration(grid[i_1,i_2], c, maxiter)
        end
    end
    mask
end

function indomain_grid(grid, js::JuliaSet)
    if isequal(left(grid),left(js.box)) && isequal(right(grid),right(js.box))
        if haskey(js.maskcache, size(grid,1))
            mask = js.maskcache[size(grid,1)]
        else # compute mask and cache it
            mask = computejuliasetgrid(grid, js.c, js.maxiter)
            js.maskcache[size(grid,1)] = mask
        end
    else # Don't cache if the grid doesn't match the bounding box
        mask = computejuliasetgrid(grid, js.c, js.maxiter)
    end
    mask
end

indomain(x::SVector{2}, js::JuliaSet) = juliasetiteration(x, js.c, js.maxiter)

show(io::IO, js::JuliaSet) = print(io, "A particular Julia Set also known as the Douady rabbit")

boundingbox(js::JuliaSet) = js.box
