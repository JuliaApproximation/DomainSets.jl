module DomainSetsGeometryBasicsExt

using DomainSetsCore, DomainSets, GeometryBasics

const GB = GeometryBasics

using StaticArrays

import DomainSets: todomainset

DomainSetsCore.DomainStyle(d::GB.AbstractGeometry) = IsDomain()
DomainSetsCore.domaineltype(d::GB.AbstractGeometry{DIM,T}) where {DIM,T} = SVector{DIM,T}


## The Point type

# repeat because AbstractPoint does not inherit from AbstractGeometry, it is a vector
DomainSetsCore.DomainStyle(d::GB.AbstractPoint) = IsDomain()
DomainSetsCore.domaineltype(d::GB.AbstractPoint{DIM,T}) where {DIM,T} = SVector{DIM,T}

DomainSets.convert_eltype(::Type{SVector{N,T}}, d::GB.Point{N,T}) where {N,T} = d
DomainSets.convert_eltype(::Type{SVector{N,T}}, d::GB.Point{N}) where {N,T} = GB.Point{N,T}(d.data)

todomainset(d::GB.AbstractPoint{N,T}) where {N,T} = DomainSets.Point{SVector{N,T}}(d)

DomainSets.canonicaldomain(::DomainSets.Equal, d::GB.AbstractPoint) = todomainset(d)

DomainSets.canonicaldomain(d::GB.AbstractPoint) = canonicaldomain(todomainset(d))
DomainSets.mapfrom_canonical(d::GB.AbstractPoint) = DomainSets.mapfrom_canonical(todomainset(d))


## The HyperRectangle primitive

DomainSets.convert_eltype(::Type{SVector{N,T}}, d::GB.HyperRectangle{N,T}) where {N,T} = d
DomainSets.convert_eltype(::Type{SVector{N,T}}, d::GB.HyperRectangle{N}) where {N,T} =
    GB.HyperRectangle{N,T}(d.origin, d.widths)

todomainset(d::GB.HyperRectangle{N,T}) where {N,T} =
    DomainSets.Rectangle{SVector{N,T}}(d.origin, d.origin+d.widths)

DomainSets.canonicaldomain(::DomainSets.Equal, d::GB.HyperRectangle) =
    todomainset(d)

DomainSets.canonicaldomain(d::GB.HyperRectangle) =
    canonicaldomain(todomainset(d))
DomainSets.mapfrom_canonical(d::GB.HyperRectangle) =
    DomainSets.mapfrom_canonical(todomainset(d))


## The HyperSphere primitive

DomainSets.convert_eltype(::Type{SVector{N,T}}, d::GB.HyperSphere{N,T}) where {N,T} = d
DomainSets.convert_eltype(::Type{SVector{N,T}}, d::GB.HyperSphere{N}) where {N,T} =
    GB.Sphere{N,T}(d.center, d.r)

todomainset(d::GB.HyperSphere{N,T}) where {N,T} =
    DomainSets.Ball{SVector{N,T}}(d.r, d.center)

DomainSets.canonicaldomain(::DomainSets.Equal, d::GB.HyperSphere) = todomainset(d)

DomainSets.canonicaldomain(d::GB.HyperSphere) = canonicaldomain(todomainset(d))
DomainSets.mapfrom_canonical(d::GB.HyperSphere) = DomainSets.mapfrom_canonical(todomainset(d))

end # module
