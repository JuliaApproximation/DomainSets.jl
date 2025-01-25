
# New deprecations in 0.8:
import Base: isreal
@deprecate isreal(d::Domain) isrealdomain(d)

LinearAlgebra.cross(d1::Domain, domains...) = linearalgebra_x_becomes_domainsets_x(d1, domains...)
@deprecate linearalgebra_x_becomes_domainsets_x(d1::Domain, domains...) cartesianproduct(d1, domains...)

####################
# Change in exports
####################

# Quite a few functions were exported in DomainSets 0.7.15 but are no longer
# exported in DomainSets 0.8. Of those, some have moved to FunctionMaps, some
# are just internal to DomainSets.

# Any missing exports not in the two lists below have been deprecated and removed.

# A) Functions that now live in FunctionMaps (exported or unexported)
# "AffineMap"
# "ConstantMap"
# "IdentityMap"
# "LinearMap"
# "Map"
# "MapRef"
# "ProductMap"
# "Translation"
# "UnityMap"
# "ZeroMap"
# "affinematrix"
# "affinevector"
# "applymap"
# "codomaintype"
# "composedmap"
# "diffvolume"
# "domaintype"
# "inverse"
# "isaffinemap"
# "isconstantmap"
# "isequalmap"
# "islinearmap"
# "jacdet"
# "jacobian"
# "leftinverse"
# "mapconstant"
# "mapsize"
# "productmap"
# "rightinverse"

# B) Functions that are no longer exported but remain available in DomainSets
# "CartToPolarMap"
# "ComposedMap"
# "PolarToCartMap"
# "StaticIdentityMap"
# "TypedMap"
# "VectorIdentityMap"
# "rotation_map"
