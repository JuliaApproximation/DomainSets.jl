# interval.jl

##################
### An interval
##################

struct Interval{T} <: Domain{1}
    a     ::  T
    b     ::  T

    Interval{T}(a = 0, b = 1) where T = new(a,b)
end

Interval() = Interval{Int}()

Interval{T}(::Type{T}) = Interval{T}()

Interval{T <: Number}(a::T, b::T) = Interval{T}(a, b)
Interval{S <: Number, T <: Number}(a::S, b::T) = Interval(promote(a,b)...)


indomain(x, d::Interval) = in(x, d.a, d.b)

left(d::Interval) = d.a
right(d::Interval) = d.b

# Arithmetic operations

(+)(d::Interval, x::Number) = Interval(d.a+x, d.b+x)

(*)(a::Number, d::Interval) = Interval(a*d.a, a*d.b)
(*)(d::Interval, a::Number) = a * d


boundingbox(d::Interval) = BBox(left(d), right(d))

show(io::IO, d::Interval) = print(io, "the interval [", d.a, ", ", d.b, "]")

const unitinterval = Interval()

function union(d1::Interval, d2::Interval)
    a = left(d1)
    b = right(d1)
    c = left(d2)
    d = right(d2)

    if (b < c) || (a > d)
        DomainUnion(d1, d2)
    else
        Interval(min(a, c), max(b, d))
    end
end

function intersect(d1::Interval, d2::Interval)
    a = left(d1)
    b = right(d1)
    c = left(d2)
    d = right(d2)

    if (b < c) || (a > d)
        EmptyDomain(Val{1}())
    else
        Interval(max(a, c), min(b, d))
    end
end
