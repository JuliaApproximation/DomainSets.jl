# characteristic.jl

"""
A domain that is described by its characteristic function.
"""
struct Characteristic{N,T} <: Domain{T}
    char    ::  Function
    box     ::  BBox{N,T}
end

Characteristic(char::Function, dom::Domain{T}) where {T} = Characteristic(char,boundingbox(dom))

indomain(x, c::Characteristic) = c.char(x)

boundingbox(c::Characteristic) = c.box

show(io::IO, c::Characteristic) = print(io, "a domain described by a characteristic function")
