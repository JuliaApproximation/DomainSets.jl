# Domains.jl

[![Build Status](https://travis-ci.org/JuliaApproximation/Domains.jl.svg?branch=master)](https://travis-ci.org/JuliaApproximation/Domains.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaApproximation/Domains.jl/badge.svg)](https://coveralls.io/github/JuliaApproximation/Domains.jl)


Domains.jl is a package designed to represent simple infinite sets, that
can be used to encode domains of functions. For example, the domain of the
function `log(x::Float64)` is the infinite open interval, which is represented
by the type `Halfline{Float64}()`.

## Examples

### Intervals

### Rectangles

### Circles

### Disks

### Union of domains

## The domain interface

A domain is any type that implements the functions `eltype` and `in`. If
`d` is an instance of a type that implements the domain interface, then
the domain consists of all `x` that is an `eltype(d)` such that `x in d`
returns true.

Domains often represent continuous mathematical domains, for example, a domain
`d`  representing the interval `[0,1]` would have `eltype(d) == Int` but still
have `0.2 in d` return true.

## The `Domain` type

Domains.jl contains an abstract type `Domain{T}`. All subtypes of `Domain{T}`
must implement the domain interface, and in addition support `convert(Domain{T}, d)`.
