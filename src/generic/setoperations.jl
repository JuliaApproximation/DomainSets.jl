# The union, intersection and difference of domains are represented with lazy domains.

Base.issubset(d1::AnyDomain, d2::AnyDomain) = issubset_domain(domain(d1), domain(d2))

issubset_domain(d1, d2) =
	promotable_domains(d1, d2) && issubset1(promote_domains(d1, d2)...)

issubset1(d1, d2) = simplifies(d1) ? issubset_domain(simplify(d1),d2) : issubset2(d1, d2)
issubset2(d1, d2) = simplifies(d2) ? issubset_domain(d1, simplify(d2)) : d1 == d2
# this last fallback is only an approximation of the truth: if d1 equals d2, then
# d1 is a subset of d2, but the reverse is not true. So we might be returning false
# even when the correct mathematical answer is true.
# What `issubset` means for the optimizations below is:
# - if true: we are sure that d1 is a subset of d2
# - if false: either it really is false, or we don't really know

issubset1(d1::Number, d2) = in(d1, d2)
issubset1(d1::AbstractArray, d2) = all(in(d2), d1)
issubset1(d1::AbstractSet, d2) = all(in(d2), d1)


############################
# The union of two domains
############################

"""
A `UnionDomain` represents the union of a set of domains.
"""
struct UnionDomain{T,DD} <: CompositeDomain{T}
	domains	::	DD
end

"""
The `UnionDomain` and `UnionDomain{T}` constructors can be invoked in three ways:
- with a list of arguments: `UnionDomain(d1, d2, ...)`
- or with an iterable list of domains: `UnionDomain(domains)`
"""
UnionDomain() = throw(ArgumentError("Can't create an empty UnionDomain."))
UnionDomain(domains...) = UnionDomain(domains)
@deprecate UnionDomain(domain::Domain) UnionDomain((domain,))
UnionDomain(domains) = _UnionDomain(promote_domains(domains))
_UnionDomain(domains) = _UnionDomain(domaineltype(first(domains)), domains)

UnionDomain{T}(domains...) where {T} = UnionDomain{T}(domains)
@deprecate UnionDomain{T}(domain::Domain) where {T} UnionDomain{T}((domain,))
UnionDomain{T}(domains) where {T} = _UnionDomain(T, convert_eltype.(T, domains))
_UnionDomain(::Type{T}, domains) where {T} = UnionDomain{T,typeof(domains)}(domains)

# The union of domains corresponds to a logical OR of their characteristic functions
composition(d::UnionDomain) = Combination()
combine(d::UnionDomain, results) = reduce(|, results)

# Make d1 ∪ d2 invoke `uniondomain` if the arguments are domains
Base.union(domains::AnyDomain...) =	uniondomain(map(domain, domains)...)

uniondomain() = emptyspace(Any)
uniondomain(d1) = d1
uniondomain(d1, d2) = uniondomain1(promote_domains(d1, d2)...)
uniondomain1(d1, d2) = simplifies(d1) ? uniondomain(simplify(d1), d2) : uniondomain2(d1, d2)
uniondomain2(d1, d2) = simplifies(d2) ? uniondomain(d1, simplify(d2)) : default_uniondomain(d1, d2)

function default_uniondomain(d1, d2)
	if issubset_domain(d2, d1)
		d1
	elseif issubset_domain(d1, d2)
		d2
	else
		UnionDomain(d1, d2)
	end
end

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
	d1 == d2 ? d1 : UnionDomain(collect(Set(components(d1)) ∪ Set(components(d2))))
uniondomain1(d1::UnionDomain, d2) = UnionDomain(components(d1)..., d2)
uniondomain2(d1, d2::UnionDomain) = UnionDomain(d1, components(d2)...)

isequaldomain(a::UnionDomain, b::UnionDomain) = Set(components(a)) == Set(components(b))
hash(d::UnionDomain, h::UInt) = hashrec("UnionDomain", Set(d.domains), h)


convert(::Type{Domain}, v::AbstractVector{<:Domain}) = UnionDomain(v)
convert(::Type{Domain}, v::AbstractSet{<:Domain}) = UnionDomain(v)
convert(::Type{Domain}, s::AbstractSet) = UnionDomain(map(Point,collect(s)))
convert(::Type{Domain{T}}, v::AbstractVector{<:Domain}) where {T} = UnionDomain{T}(v)
convert(::Type{Domain{T}}, v::AbstractSet{<:Domain}) where {T} = UnionDomain{T}(v)
convert(::Type{Domain{T}}, s::AbstractSet) where {T} = UnionDomain{T}(map(Point,collect(s)))

similardomain(d::UnionDomain, ::Type{T}) where {T} =
    UnionDomain(convert_eltype.(T, components(d)))



choice(d::UnionDomain) = convert(eltype(d), choice(component(d,1)))

isempty(d::UnionDomain) = all(isempty, d.domains)

interior(d::UnionDomain) = UnionDomain(map(interior, components(d)))
closure(d::UnionDomain) = UnionDomain(map(closure, components(d)))

boundingbox(d::UnionDomain) = unionbox(map(boundingbox, components(d))...)

boundary(d::UnionDomain) = uniondomain(map(boundary, components(d))...)

Display.combinationsymbol(d::UnionDomain) = Display.Symbol('∪')
Display.displaystencil(d::UnionDomain) = composite_displaystencil(d)
show(io::IO, mime::MIME"text/plain", d::UnionDomain) = Display.composite_show(io, mime, d)
show(io::IO, d::UnionDomain) = Display.composite_show_compact(io, d)

## algebra

# preserve the uniondomain when mapping
map_domain(map, domain::UnionDomain) = UnionDomain(map_domain.(Ref(map), components(domain)))
mapped_domain(map, domain::UnionDomain) = UnionDomain(mapped_domain.(Ref(map), components(domain)))


# TODO: what is the correct semantics for these functions? Should we have them?
Base.minimum(d::UnionDomain) = mapreduce(minimum, min, components(d))
Base.maximum(d::UnionDomain) = mapreduce(maximum, max, components(d))
infimum(d::UnionDomain) = mapreduce(infimum, min, components(d))
supremum(d::UnionDomain) = mapreduce(supremum, max, components(d))

setdiffdomain(d1::UnionDomain, d2::UnionDomain) =
	UnionDomain(setdiffdomain.(components(d1), Ref(d2)))
function setdiffdomain1(d1::UnionDomain, d2)
	if d2 ∈ components(d1)
		el = filter(x->x!=d2, components(d1))
		if length(el) == 0
			EmptySpace{domaineltype(d1)}()
		elseif length(el) == 1
			first(el)
		else
			UnionDomain(el)
		end
	else
		UnionDomain(setdiffdomain.(components(d1), Ref(d2)))
	end
end
function setdiffdomain2(d1, d2::UnionDomain)
    result = d1
    for d in components(d2)
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
struct IntersectDomain{T,DD} <: CompositeDomain{T}
    domains ::  DD
end

"""
The `IntersectDomain` constructor can be invoked in one of three ways:
- with a list of arguments: IntersectDomain(d1, d2, ...)
- or with any iterable list of domains: IntersectDomain(domains)
"""
IntersectDomain(domains...) = IntersectDomain(domains)
@deprecate IntersectDomain(domain::Domain) IntersectDomain((domain,))
IntersectDomain(domains) = _IntersectDomain(promote_domains(domains))
_IntersectDomain(domains) = IntersectDomain{domaineltype(first(domains))}(domains)

IntersectDomain{T}(domains...) where {T} = IntersectDomain{T}(domains)
@deprecate IntersectDomain{T}(domain::Domain) where T IntersectDomain((domain,))
IntersectDomain{T}(domains) where {T} = _IntersectDomain(T, convert_eltype.(T, domains))
_IntersectDomain(::Type{T}, domains) where {T} = IntersectDomain{T,typeof(domains)}(domains)

# The intersection of domains corresponds to a logical AND of their characteristic functions
composition(d::IntersectDomain) = Combination()
combine(d::IntersectDomain, results) = reduce(&, results)


# Make d1 ∩ d2 invoke `intersectdomain` if the arguments are Domains
Base.intersect(domains::AnyDomain...) = intersectdomain(map(domain, domains)...)

intersectdomain() = emptyspace(Any)
intersectdomain(d1) = d1
intersectdomain(d1, d2) = intersectdomain1(promote_domains(d1, d2)...)

intersectdomain1(d1, d2) = simplifies(d1) ? intersectdomain(simplify(d1), d2) : intersectdomain2(d1, d2)
intersectdomain2(d1, d2) = simplifies(d2) ? intersectdomain(d1, simplify(d2)) : default_intersectdomain(d1, d2)

function default_intersectdomain(d1, d2)
	if issubset_domain(d2, d1)
		d2
	elseif issubset_domain(d1, d2)
		d1
	else
		IntersectDomain(d1, d2)
	end
end

intersectdomain(d1, d2, d3) = _id3(promote_domains(d1, d2, d3)...)
_id3(d1, d2, d3) =
	_id3(d1, d2, d3, intersectdomain(d1, d2), intersectdomain(d2, d3), intersectdomain(d1, d3))
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
	d1 == d2 ? d1 : IntersectDomain(components(d1)..., components(d2)...)
intersectdomain1(d1::IntersectDomain, d2) = IntersectDomain(components(d1)..., d2)
intersectdomain2(d1, d2::IntersectDomain) = IntersectDomain(d1, components(d2)...)


function intersectdomain(d1::UnionDomain, d2::UnionDomain)
    d1 == d2 && return d1
    uniondomain(intersectdomain.(Ref(d1), components(d2))...)
end
intersectdomain1(d1::UnionDomain, d2) = uniondomain(intersectdomain.(d1.domains, Ref(d2))...)
intersectdomain2(d1, d2::UnionDomain) = uniondomain(intersectdomain.(Ref(d1), d2.domains)...)

Base.:&(d1::AnyDomain, d2::AnyDomain) = intersectdomain(domain(d1),domain(d2))

function intersectdomain(d1::ProductDomain, d2::ProductDomain)
	if compatibleproductdims(d1, d2)
        ProductDomain(map(intersectdomain, components(d1), components(d2)))
    else
        IntersectDomain(d1, d2)
    end
end

similardomain(d::IntersectDomain, ::Type{T}) where {T} =
    IntersectDomain(convert_eltype.(T, components(d)))

isequaldomain(a::IntersectDomain, b::IntersectDomain) = Set(components(a)) == Set(components(b))
hash(d::IntersectDomain, h::UInt) = hashrec("IntersectDomain", Set(components(d)), h)

boundingbox(d::IntersectDomain) = intersectbox(map(boundingbox, components(d))...)

Display.combinationsymbol(d::IntersectDomain) = Display.Symbol('∩')
Display.displaystencil(d::IntersectDomain) = composite_displaystencil(d)
show(io::IO, mime::MIME"text/plain", d::IntersectDomain) = Display.composite_show(io, mime, d)
show(io::IO, d::IntersectDomain) = Display.composite_show_compact(io, d)


#########################################
### The difference between two domains
#########################################


"A `SetdiffDomain` represents the difference between two domains."
struct SetdiffDomain{T,DD} <: CompositeDomain{T}
    domains	::	DD
	function SetdiffDomain{T,DD}(domains::DD) where {T,DD}
		@assert length(domains) == 2
		new(domains)
	end
end

SetdiffDomain(d1, d2) = _SetdiffDomain(promote_domains((d1, d2))...)
_SetdiffDomain(d1, d2) = SetdiffDomain{domaineltype(d1)}((d1,d2))
SetdiffDomain{T}(domains) where {T} = SetdiffDomain{T,typeof(domains)}(domains)

# The difference between two domains corresponds to a logical AND NOT of their characteristic functions
composition(d::SetdiffDomain) = Combination()
combine(d::SetdiffDomain, results) = results[1] & !results[2]

# It is difficult to calculate approximate membership exactly, but we can at
# least not enlarge the subtracted domain by invoking in rather than approx_in on it.
_approx_indomain(x, d::SetdiffDomain, comp::Combination, domains, tolerance) =
    approx_in(x, domains[1], tolerance) & !in(x, domains[2])

similardomain(d::SetdiffDomain, ::Type{T}) where {T} =
    SetdiffDomain(convert_eltype(T, d.domains[1]), convert_eltype(T, d.domains[2]))

# use \ as a synomym for setdiff, in the context of domains (though, generically,
# \ means left division in Julia)
Base.:\(d1::AnyDomain, d2::AnyDomain) = setdiffdomain(domain(d1), domain(d2))

# Make setdiff invoke `setdiffdomain` if the arguments are domains
setdiff(d1::AnyDomain, d2::AnyDomain) = setdiffdomain(domain(d1), domain(d2))

setdiffdomain(d1, d2) = setdiffdomain1(promote_domains(d1, d2)...)
setdiffdomain1(d1, d2) = simplifies(d1) ? setdiffdomain(simplify(d1), d2) : setdiffdomain2(d1, d2)
setdiffdomain2(d1, d2) = simplifies(d2) ? setdiffdomain(d1, simplify(d2)) : default_setdiffdomain(d1, d2)

function default_setdiffdomain(d1, d2)
	if isempty(d2)
		d1
	elseif issubset_domain(d1,d2)
		emptyspace(d1)
	else
		SetdiffDomain(d1, d2)
	end
end

# avoid nested difference domains
setdiffdomain1(d1::SetdiffDomain, d2) = setdiffdomain(d1.domains[1], uniondomain(d2, d1.domains[2]))

isequaldomain(a::SetdiffDomain, b::SetdiffDomain) = a.domains == b.domains

boundingbox(d::SetdiffDomain) =  boundingbox(d.domains[1])

Display.combinationsymbol(d::SetdiffDomain) = Display.Symbol('\\')
Display.displaystencil(d::SetdiffDomain) = composite_displaystencil(d)
show(io::IO, mime::MIME"text/plain", d::SetdiffDomain) = Display.composite_show(io, mime, d)
show(io::IO, d::SetdiffDomain) = Display.composite_show_compact(io, d)
