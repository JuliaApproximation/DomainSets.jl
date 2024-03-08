using FunctionMaps

# for test_common.jl
using FunctionMaps:
    convert_numtype,
    convert_prectype,
    promote_numtype,
    promote_prectype,
    euclideandimension

# for test_interface.jl
using FunctionMaps:
    MapStyle, IsMap, NotMap,
    functionmap, checkmap

# for test_generic.jl
using FunctionMaps:
    convert_domaintype, convert_codomaintype,
    map_hash,
    LazyInverse, jacobian!

# for test_maps.jl
using FunctionMaps: ScalarAffineMap,
    VectorAffineMap,
    StaticAffineMap,
    GenericAffineMap,
    ScalarLinearMap,
    VectorLinearMap,
    StaticLinearMap,
    GenericLinearMap,
    ScalarTranslation,
    VectorTranslation,
    StaticTranslation,
    GenericTranslation,
    ProductMap,
    TupleProductMap, VcatMap, VectorProductMap,
    WrappedMap,
    interval_map, multiply_map,
    SumMap, sum_map,
    composedmap, ComposedMap,
    composite_jacobian, sum_jacobian,
    CartToPolarMap, PolarToCartMap,
    UnitCircleMap, AngleMap, UnitDiskMap,
    VectorToComplex

# for test_affine.jl
using FunctionMaps:
    to_matrix, to_vector,
    matrix_pinv

# for test_arithmetics.jl
using FunctionMaps:
    interval_map, bounded_interval_map
