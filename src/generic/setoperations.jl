# The union, intersection and difference of domains are represented with lazy domains.

issubset(d1::Domain, d2::Domain) = compatible_eltype(d1, d2) && issubset1(promote_domains(d1, d2)...)
issubset(d1::Domain, d2) = issubset1(promote_domains(d1, d2)...)
issubset(d1, d2::Domain) = issubset1(promote_domains(d1, d2)...)
issubset1(d1, d2) = issubset2(d1, d2)
issubset2(d1, d2) = d1 == d2


############################
# The union of two domains
############################

"""
```using DomainSets```

A `UnionDomain` represents the union of a set of domains.
"""
struct UnionDomain{T,DD} <: CompositeLazyDomain{T}
	domains	::	DD
end

"""
The `UnionDomain` and `UnionDomain{T}` constructors can be invoked in three ways:
- with a list of arguments: UnionDomain(d1, d2, ...)
- with a single domain: UnionDomain(d::Domain)
- or with any iterable list of domains: UnionDomain(domains)
"""
UnionDomain(domains...) = UnionDomain(domains)
UnionDomain(d::Domain) = UnionDomain((d,))
UnionDomain(domains) = _UnionDomain(promote_domains(domains))
_UnionDomain(domains) = _UnionDomain(eltype(first(domains)), domains)

UnionDomain{T}(domains...) where {T} = UnionDomain{T}(domains)
UnionDomain{T}(d::Domain) where {T} = UnionDomain{T}((d,))
UnionDomain{T}(domains) where {T} = _UnionDomain(T, convert_eltype.(T, domains))
_UnionDomain(::Type{T}, domains) where {T} = UnionDomain{T,typeof(domains)}(domains)

# The union of domains corresponds to a logical OR of their characteristic functions
composition(d::UnionDomain) = Combination()
combine(d::UnionDomain, results) = reduce(|, results)

# Make d1 ∪ d2 invoke `uniondomain` if one of the first two arguments is a Domain
union(d1::Domain, d2::Domain, domains...) = uniondomain(d1, d2, domains...)
union(d1::Domain, d2, domains...) = uniondomain(d1, d2, domains...)
union(d1, d2::Domain, domains...) = uniondomain(d1, d2, domains...)

uniondomain() = EmptySpace{Any}()
uniondomain(d1) = d1
uniondomain(d1, d2) = uniondomain1(promote_domains(d1, d2)...)
uniondomain1(d1, d2) = issubset(d1, d2) ? d2 : uniondomain2(d1, d2)
uniondomain2(d1, d2) = issubset(d2, d1) ? d1 : UnionDomain(d1, d2)

uniondomain(d1, d2, d3) = _ud3(promote_domains(d1, d2, d3)...)
_ud3(d1, d2, d3) =
	_ud3(d1, d2, d3, uniondomain(d1, d2), uniondomain(d2, d3), uniondomain(d1, d3))
_ud3(d1, d2, d3, d12::UnionDomain, d23::UnionDomain, d13::UnionDomain) =
	uniondomain(d12, d3)
_ud3(d1, d2, d3, d12, d23::UnionDomain, d13::UnionDomain) =
	uniondomain(d12, d3)
_ud3(d1, d2, d3, d12::UnionDomain, d23, d13::UnionDomain) =
	uniondomain(d23, d1)
_ud3(d1, d2, d3, d12::UnionDomain, d23::UnionDomain, d13) =
	uniondomain(d13, d2)
_ud3(d1, d2, d3, d12::UnionDomain, d23, d13) =
	uniondomain(d23, d1)
_ud3(d1, d2, d3, d12, d23::UnionDomain, d13) =
	uniondomain(d12, d3)
_ud3(d1, d2, d3, d12, d23, d13::UnionDomain) =
	uniondomain(d12, d3)
_ud3(d1, d2, d3, d12, d23, d13) =
	uniondomain(d12, d3)

uniondomain(d1, d2, d3, domains...) = _ud(promote_domains(d1, d2, d3, domains...)...)
_ud(d1, d2, d3, domains...) = uniondomain(uniondomain(d1, d2, d3), domains...)

# avoid nested union domains
uniondomain(d1::UnionDomain, d2::UnionDomain) =
	d1 == d2 ? d1 : UnionDomain(elements(d1)..., elements(d2)...)
uniondomain1(d1::UnionDomain, d2) = UnionDomain(elements(d1)..., d2)
uniondomain2(d1, d2::UnionDomain) = UnionDomain(d1, elements(d2)...)

==(a::UnionDomain, b::UnionDomain) = Set(elements(a)) == Set(elements(b))


convert(::Type{Domain}, v::AbstractVector{<:Domain}) = UnionDomain(v)
convert(::Type{Domain}, v::AbstractSet{<:Domain}) = UnionDomain(v)
convert(::Type{Domain}, s::AbstractSet) = UnionDomain(map(Point,collect(s)))
convert(::Type{Domain{T}}, v::AbstractVector{<:Domain}) where {T} = UnionDomain{T}(v)
convert(::Type{Domain{T}}, v::AbstractSet{<:Domain}) where {T} = UnionDomain{T}(v)
convert(::Type{Domain{T}}, s::AbstractSet) where {T} = UnionDomain{T}(map(Point,collect(s)))

similardomain(d::UnionDomain, ::Type{T}) where {T} =
    UnionDomain(convert.(Domain{T}, elements(d)))

hash(d::UnionDomain, h::UInt) = hash(Set(d.domains), h)


point_in_domain(d::UnionDomain) = convert(eltype(d), point_in_domain(element(d,1)))

isempty(d::UnionDomain) = all(isempty, d.domains)

interior(d::UnionDomain) = UnionDomain(map(interior, elements(d)))
closure(d::UnionDomain) = UnionDomain(map(closure, elements(d)))

function show(io::IO, d::UnionDomain)
    print(io, "the union of $(numelements(d)) domains:\n")
    for (i,e) in enumerate(elements(d))
        print(io, "\t$i.\t: ", e, "\n")
    end
end


## ualgebra

# preserve the uniondomain when mapping
map_domain(map, domain::UnionDomain) = UnionDomain(map_domain.(Ref(map), elements(domain)))
mapped_domain(map, domain::UnionDomain) = UnionDomain(mapped_domain.(Ref(map), elements(domain)))

# for op in (:+, :-, :*, :/)
#     @eval begin
#         $op(domain::UnionDomain, x::Number) = UnionDomain(broadcast($op, elements(domain), x))
#         $op(x::Number, domain::UnionDomain) = UnionDomain(broadcast($op, x, elements(domain)))
#     end
# end
#
# \(x::Number, domain::UnionDomain) = UnionDomain(broadcast(\, x, elements(domain)))


for (op, mop) in ((:minimum, :min), (:maximum, :max), (:infimum, :min), (:supremum, :max))
    @eval $op(d::UnionDomain) = mapreduce($op, $mop, elements(d))
end


setdiffdomain(d1::UnionDomain, d2::UnionDomain) =
	UnionDomain(setdiffdomain.(elements(d1), Ref(d2)))

function setdiffdomain(d1::UnionDomain, d2::Domain)
    s = Set(elements(d1))
    # check if any element is in d1 and just remove
    s2 = Set(setdiff(s, tuple(d2)))
    s2 ≠ s && return UnionDomain(s2)
    UnionDomain(setdiffdomain.(elements(d1), Ref(d2)))
end

function setdiffdomain(d1::Domain, d2::UnionDomain)
    result = d1
    for d in elements(d2)
        result = setdiffdomain(result, d)
    end
    result
end


##############################
# The intersection of domains
##############################


"""
An `IntersectDomain` represents the intersection of a set of domains.
"""
struct IntersectDomain{T,DD} <: CompositeLazyDomain{T}
    domains ::  DD
end

"""
The `IntersectDomain` constructor can be invoked in one of three ways:
- with a list of arguments: IntersectDomain(d1, d2, ...)
- with a single domain: IntersectDomain(d::Domain)
- or with any iterable list of domains: IntersectDomain(domains)
"""
IntersectDomain(domains...) = IntersectDomain(domains)
IntersectDomain(d::Domain) = IntersectDomain((d,))
IntersectDomain(domains) = _IntersectDomain(promote_domains(domains))
_IntersectDomain(domains) = IntersectDomain{eltype(first(domains))}(domains)

IntersectDomain{T}(domains...) where {T} = IntersectDomain{T}(domains)
IntersectDomain{T}(d::Domain) where {T} = IntersectDomain{T}((d,))
IntersectDomain{T}(domains) where {T} = _IntersectDomain(T, convert_eltype.(T, domains))
_IntersectDomain(::Type{T}, domains) where {T} = IntersectDomain{T,typeof(domains)}(domains)

# The intersection of domains corresponds to a logical AND of their characteristic functions
composition(d::IntersectDomain) = Combination()
combine(d::IntersectDomain, results) = reduce(&, results)


# Make d1 ∩ d2 invoke `intersectdomain` if one of the first two arguments is a Domain
intersect(d1::Domain, d2::Domain, domains...) = intersectdomain(d1, d2, domains...)
intersect(d1::Domain, d2, domains...) = intersectdomain(d1, d2, domains...)
intersect(d1, d2::Domain, domains...) = intersectdomain(d1, d2, domains...)

intersectdomain() = EmptySpace{Any}()
intersectdomain(d1) = d1
intersectdomain(d1, d2) = intersectdomain1(promote_domains(d1, d2)...)
intersectdomain1(d1, d2) = issubset(d1, d2) ? d1 : intersectdomain2(d1, d2)
intersectdomain2(d1, d2) = issubset(d2, d1) ? d2 : IntersectDomain(d1, d2)

intersectdomain(d1, d2, d3) = _id3(promote_domains(d1, d2, d3)...)
_id3(d1, d2, d3) =
	_id3(d1, d2, d3, intersectdomain(d1, d2), intersectdomain(d2, d3), uniondomain(d1, d3))
_id3(d1, d2, d3, d12::IntersectDomain, d23::IntersectDomain, d13::IntersectDomain) =
	intersectdomain(d12, d3)
_id3(d1, d2, d3, d12, d23::IntersectDomain, d13::IntersectDomain) =
	intersectdomain(d12, d3)
_id3(d1, d2, d3, d12::IntersectDomain, d23, d13::IntersectDomain) =
	intersectdomain(d23, d1)
_id3(d1, d2, d3, d12::IntersectDomain, d23::IntersectDomain, d13) =
	intersectdomain(d13, d2)
_id3(d1, d2, d3, d12::IntersectDomain, d23, d13) =
	intersectdomain(d23, d1)
_id3(d1, d2, d3, d12, d23::IntersectDomain, d13) =
	intersectdomain(d12, d3)
_id3(d1, d2, d3, d12, d23, d13::IntersectDomain) =
	intersectdomain(d12, d3)
_id3(d1, d2, d3, d12, d23, d13) =
	intersectdomain(d12, d3)

intersectdomain(d1, d2, d3, domains...) = _id(promote_domains(d1, d2, d3, domains...)...)
_id(d1, d2, d3, domains...) = intersectdomain(intersectdomain(d1, d2, d3), domains...)

# avoid nested intersect domains
intersectdomain(d1::IntersectDomain, d2::IntersectDomain) =
	d1 == d2 ? d1 : IntersectDomain(elements(d1)..., elements(d2)...)
intersectdomain1(d1::IntersectDomain, d2) = IntersectDomain(elements(d1)..., d2)
intersectdomain2(d1, d2::IntersectDomain) = IntersectDomain(d1, elements(d2)...)


function intersectdomain(d1::UnionDomain, d2::UnionDomain)
    d1 == d2 && return d1
    uniondomain(intersectdomain.(Ref(d1), elements(d2))...)
end
intersectdomain(d1::UnionDomain, d2::Domain) = uniondomain(intersectdomain.(d1.domains, Ref(d2))...)
intersectdomain(d1::Domain, d2::UnionDomain) = uniondomain(intersectdomain.(Ref(d1), d2.domains)...)

(&)(d1::Domain, d2::Domain) = intersectdomain(d1,d2)

function intersectdomain(d1::ProductDomain, d2::ProductDomain)
	if compatibleproduct(d1, d2)
        ProductDomain(map(intersectdomain, elements(d1), elements(d2)))
    else
        IntersectDomain(d1, d2)
    end
end

similardomain(d::IntersectDomain, ::Type{T}) where {T} =
    IntersectDomain(convert.(Domain{T}, elements(d)))

==(a::IntersectDomain, b::IntersectDomain) = Set(elements(a)) == Set(elements(b))

function show(io::IO, d::IntersectDomain)
    print(io, "the intersection of $(numelements(d)) domains:\n")
	for i in 1:numelements(d)
    	print(io, "\t$(i).\t: ", element(d,i), "\n")
	end
end


#########################################
### The difference between two domains
#########################################

@deprecate DifferenceDomain SetdiffDomain

"A `SetdiffDomain` represents the difference between two domains."
struct SetdiffDomain{T,DD} <: CompositeLazyDomain{T}
    domains	::	DD
	function SetdiffDomain{T,DD}(domains::DD) where {T,DD}
		@assert length(domains) == 2
		new(domains)
	end
end

SetdiffDomain(d1, d2) = _SetdiffDomain(promote_domains((d1, d2))...)
_SetdiffDomain(d1, d2) = SetdiffDomain{eltype(d1)}((d1,d2))
SetdiffDomain{T}(domains) where {T} = SetdiffDomain{T,typeof(domains)}(domains)

# The difference between two domains corresponds to a logical AND NOT of their characteristic functions
composition(d::SetdiffDomain) = Combination()
combine(d::SetdiffDomain, results) = results[1] & !results[2]

# It is difficult to calculate approximate membership exactly, but we can at
# least not enlarge the subtracted domain by invoking in rather than approx_in on it.
_approx_indomain(x, d::SetdiffDomain, comp::Combination, domains, tolerance) =
    approx_in(x, domains[1], tolerance) & !in(x, domains[2])

similardomain(d::SetdiffDomain, ::Type{T}) where {T} =
    SetdiffDomain(convert(Domain{T}, d.domains[1]), convert(Domain{T}, d.domains[2]))

# use \ as a synomym for setdiff, in the context of domains (though, generically,
# \ means left division in Julia)
\(d1::Domain, d2::Domain) = setdiffdomain(d1, d2)
\(d1::Domain, d2) = setdiffdomain(d1, d2)
\(d1, d2::Domain) = setdiffdomain(d1, d2)

# Make setdiff invoke `setdiffdomain` if one of the arguments is a Domain
setdiff(d1::Domain, d2::Domain) = setdiffdomain(d1, d2)
setdiff(d1::Domain, d2) = setdiffdomain(d1, d2)
setdiff(d1, d2::Domain) = setdiffdomain(d1, d2)

setdiffdomain(d1, d2) = setdiffdomain1(promote_domains(d1, d2)...)
setdiffdomain1(d1, d2) = setdiffdomain2(d1, d2)
function setdiffdomain2(d1, d2)
	if isempty(d2)
		d1
	elseif issubset(d1,d2)
		EmptySpace{eltype(d1)}()
	else
		SetdiffDomain(d1, d2)
	end
end

# avoid nested difference domains
setdiffdomain1(d1::SetdiffDomain, d2) = setdiffdomain(d1.domains[1], uniondomain(d2, d1.domains[2]))

==(a::SetdiffDomain, b::SetdiffDomain) = a.domains == b.domains

function show(io::IO, d::SetdiffDomain)
    print(io, "the difference of 2 domains:\n")
    print(io, "\t1.\t: ", element(d, 1), "\n")
    print(io, "\t2.\t: ", element(d, 2), "\n")
end
