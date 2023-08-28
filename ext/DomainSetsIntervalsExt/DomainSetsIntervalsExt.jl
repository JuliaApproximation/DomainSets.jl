module DomainSetsIntervalsExt

using DomainSetsCore, DomainSets, Intervals

import DomainSets: todomainset

## The Intervals.Interval type

DomainSetsCore.DomainStyle(d::Intervals.AbstractInterval) = IsDomain()
DomainSetsCore.domaineltype(d::Intervals.Interval) = eltype(d)

DomainSets.convert_eltype(::Type{T}, d::Intervals.Interval{T}) where {T} = d
DomainSets.convert_eltype(::Type{T}, d::Intervals.Interval) where {T} =
    Intervals.Interval(T(d.first), T(d.last))

todomainset(d::Intervals.Interval{T,Closed,Closed}) where T =
    DomainSets.Interval{:closed,:closed,T}(d.first, d.last)
todomainset(d::Intervals.Interval{T,Closed,Open}) where T =
    DomainSets.Interval{:closed,:open,T}(d.first, d.last)
todomainset(d::Intervals.Interval{T,Open,Closed}) where T =
    DomainSets.Interval{:open,:closed,T}(d.first, d.last)
todomainset(d::Intervals.Interval{T,Open,Open}) where T =
    DomainSets.Interval{:open,:open,T}(d.first, d.last)

DomainSets.canonicaldomain(::DomainSets.Equal, d::Intervals.Interval) =
    todomainset(d)

DomainSets.canonicaldomain(d::Intervals.Interval) =
    canonicaldomain(todomainset(d))
DomainSets.mapfrom_canonical(d::Intervals.Interval) =
    DomainSets.mapfrom_canonical(todomainset(d))


## The Intervals.IntervalSet type

DomainSetsCore.DomainStyle(d::Intervals.IntervalSet) = IsDomain()
DomainSetsCore.domaineltype(d::Intervals.IntervalSet) = domaineltype(d.items[1])

DomainSets.convert_eltype(::Type{T}, d::Intervals.IntervalSet) where T =
    Intervals.IntervalSet(map(d->DomainSets.convert_eltype(T, d), d.items))

todomainset(d::Intervals.IntervalSet) = UnionDomain(map(todomainset, d.items))

DomainSets.canonicaldomain(::DomainSets.Equal, d::Intervals.IntervalSet) =
    todomainset(d)

DomainSets.canonicaldomain(d::Intervals.IntervalSet) =
    canonicaldomain(todomainset(d))
DomainSets.mapfrom_canonical(d::Intervals.IntervalSet) =
    DomainSets.mapfrom_canonical(todomainset(d))

end # module
