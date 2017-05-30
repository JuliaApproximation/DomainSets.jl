# arithmetics.jl
# Code for computing with domains.

# Make sure domains only need to implement addition/multiplication with numbers to the right
(+)(x::Number, d::Domain) = d + x
(+)(x::AbstractVector, d::Domain) = d + x
(*)(x::Number, d::Domain) = d * x

(/)(d::Domain, x::Number) = d * (1/x)


################################################################################
### The union of two domains
################################################################################

struct DomainUnion{D1,D2,N} <: Domain{N}
    d1    ::  D1
    d2    ::  D2

    DomainUnion{D1,D2,N}(d1::Domain{N}, d2::Domain{N}) where {D1,D2,N} = new(d1, d2)
end

DomainUnion{N}(d1::Domain{N}, d2::Domain{N}) = DomainUnion{typeof(d1),typeof(d2),N}(d1, d2)

union(d1::Domain, d2::Domain) = (d1 == d2 ? d1 : DomainUnion(d1,d2))


# The union of two domains corresponds to a logical OR of their characteristic functions
indomain(x, d::DomainUnion) = in(x, d.d1) || in(x, d.d2)

function indomain_grid(grid, d::DomainUnion)
    z1 = indomain_grid(grid, d.d1)
    z2 = indomain_grid(grid, d.d2)
    z1 .| z2
end

(+)(d1::Domain, d2::Domain) = union(d1,d2)
(|)(d1::Domain, d2::Domain) = union(d1,d2)


boundingbox(d::DomainUnion) = boundingbox(d.d1) ∪ boundingbox(d.d2)

function show(io::IO, d::DomainUnion)
    print(io, "a union of two domains: \n")
    print(io, "    First domain: ", d.d1, "\n")
    print(io, "    Second domain: ", d.d2, "\n")
end


################################################################################
### The intersection of two domains
################################################################################

struct DomainIntersection{D1,D2,N} <: Domain{N}
    d1    ::  D1
    d2    ::  D2

    DomainIntersection{D1,D2,N}(d1::Domain{N}, d2::Domain{N}) where {D1,D2,N} = new(d1, d2)
end

DomainIntersection{N}(d1::Domain{N},d2::Domain{N}) = DomainIntersection{typeof(d1),typeof(d2),N}(d1, d2)

# The intersection of two domains corresponds to a logical AND of their characteristic functions
indomain(x, d::DomainIntersection) = in(x, d.d1) && in(x, d.d2)

function indomain_grid(grid, d::DomainIntersection)
    z1 = indomain_grid(grid, d.d1)
    z2 = indomain_grid(grid, d.d2)
    z1 .& z2
end

(&)(d1::Domain, d2::Domain) = intersect(d1,d2)

intersect(d1::Domain, d2::Domain) = (d1 == d2 ? d1 : DomainIntersection(d1,d2))

function intersect(d1::ProductDomain, d2::ProductDomain)
    @assert ndims(d1) == ndims(d2)
    if nb_elements(d1) == nb_elements(d2)
        Product([intersect(element(d1,i), element(d2,i)) for i in 1:nb_elements(d1)]...)
    else
        DomainIntersection(d1, d2)
    end
end



boundingbox(d::DomainIntersection) = boundingbox(d.d1) ∩ boundingbox(d.d2)

function show(io::IO, d::DomainIntersection)
    print(io, "the intersection of two domains: \n")
    print(io, "    First domain: ", d.d1, "\n")
    print(io, "    Second domain: ", d.d1, "\n")
end


################################################################################
### The difference between two domains
################################################################################

struct DomainDifference{D1,D2,N} <: Domain{N}
    d1    ::  D1
    d2    ::  D2

    DomainDifference{D1,D2,N}(d1::Domain{N}, d2::Domain{N}) where {D1,D2,N} = new(d1, d2)
end

DomainDifference{N}(d1::Domain{N}, d2::Domain{N}) = DomainDifference{typeof(d1),typeof(d2),N}(d1,d2)

setdiff(d1::Domain, d2::Domain) = DomainDifference(d1, d2)

# The difference between two domains corresponds to a logical AND NOT of their characteristic functions
indomain(x, d::DomainDifference) = indomain(x, d.d1) && (~indomain(x, d.d2))

function indomain_grid(grid, d::DomainDifference)
    z1 = indomain_grid(grid, d.d1)
    z2 = indomain_grid(grid, d.d2)
    z1 .& (~z2)
end

(-)(d1::Domain, d2::Domain) = setdiff(d1, d2)
(\ )(d1::Domain, d2::Domain) = setdiff(d1, d2)


boundingbox(d::DomainDifference) = boundingbox(d.d1)

function show(io::IO, d::DomainDifference)
    print(io, "the difference of two domains: \n")
    print(io, "    First domain: ", d.d1, "\n")
    print(io, "    Second domain: ", d.d2, "\n")
end


################################################################################
### A revolved domain is a 2D-domain rotated about the X-axis
################################################################################

struct RevolvedDomain{D} <: Domain{3}
    d     ::  D
end

revolve(d::Domain{2}) = RevolvedDomain(d)

function indomain(x, d::RevolvedDomain)
    r = sqrt(x[2]^2+x[3])
    phi = atan2(x[2]/x[1])
    theta = acos(x[3]/r)
    indomain((x[1],r), d.d)
end


boundingbox(d::RevolvedDomain) = BBox((left(d.d)[1],left(d.d)...),(right(d.d)[1],right(d.d)...))

function show(io::IO, r::RevolvedDomain)
    print(io, "the revolution of: ", r.d1)
end


################################################################################
### A rotated domain
################################################################################

struct RotatedDomain{D,T,N,L} <: Domain{N}
   d                 ::  D
   angle             ::  Vector{T}
   rotationmatrix    ::  SMatrix{N,N,T,L}

    RotatedDomain{D,T,N,L}(d,angle,rotationmatrix) where {D,T,N,L} = new(d, angle, rotationmatrix)
end

RotatedDomain{N,T}(d::Domain{N}, angle::Vector{T}, m::SMatrix{N,N,T} = rotationmatrix(angle)) =
   RotatedDomain{typeof(d),T,N,N*N}(d, angle, m)

RotatedDomain(d::Domain{2}, theta::Number) = RotatedDomain{typeof(d),typeof(theta),2,4}(d, [theta], rotationmatrix(theta))
# types annotated to remove ambiguity
RotatedDomain{T,D}(d::D, phi::T, theta::T, psi::T) = RotatedDomain{3,T,D}(d, [phi,theta,psi], rotationmatrix(phi,theta,psi))

rotate{T}(d::Domain{2}, theta::T) = RotatedDomain(d, theta)
rotate{T}(d::Domain{3}, phi::T, theta::T, psi::T) = RotatedDomain(d, phi, theta, psi)

indomain(x, d::RotatedDomain) = indomain(d.rotationmatrix*x, d.d)

(==)(d1::RotatedDomain, d2::RotatedDomain) = (d1.d == d2.d) && (d1.angle == d2.angle) #&& (d1.rotationmatrix == d2.rotationmatrix)

# very crude bounding box (doesn't work!!!)
boundingbox(r::RotatedDomain)= sqrt(2)*boundingbox(r.d)




###########################
### A translated domain
###########################

struct TranslatedDomain{D,T,N} <: Domain{N}
    domain  ::  D
    trans   ::  SVector{N,T}

    TranslatedDomain{D,T,N}(domain::Domain{N}, trans) where {D,T,N} = new(domain, trans)
end

TranslatedDomain{N}(domain::Domain{N}, trans::SVector{N}) = TranslatedDomain{typeof(domain),eltype(trans),N}(domain, trans)

domain(d::TranslatedDomain) = d.domain

translationvector(d::TranslatedDomain) = d.trans

function indomain(x, d::TranslatedDomain)
    indomain(x-d.trans, d.domain)
end

(+)(d::Domain, trans::SVector) = TranslatedDomain(d, trans)
(+)(d::TranslatedDomain, trans::SVector) = TranslatedDomain(domain(d), trans+translationvector(d))

boundingbox(d::TranslatedDomain) = boundingbox(domain(d)) + translationvector(d)
