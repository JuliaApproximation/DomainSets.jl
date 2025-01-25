
import Base: isreal
@deprecate isreal(d::Domain) isrealdomain(d)

LinearAlgebra.cross(d1::Domain, domains...) = linearalgebra_x_becomes_domainsets_x(d1, domains...)
@deprecate linearalgebra_x_becomes_domainsets_x(d1::Domain, domains...) cartesianproduct(d1, domains...)
