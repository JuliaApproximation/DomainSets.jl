
# In this file we list all functions defined in FunctionMaps that
# were defined in the namespace of DomainSets

using FunctionMaps:
	MapRef,
	# canonical.jl
	CanonicalType, CanonicalExtensionType,
	canonicalextensiontype,
	Equal, Equivalent,
	# composite.jl
	ComposedMap, similarmap, codomaintype, applymap,
	applymap_rec, mapsize, jacobian, backpropagate,
	inverse, leftinverse, rightinverse,
	leftinverse_rec, rightinverse_rec,
	composedmap,
	isequalmap, map_hash,
	MulMap, multiply_map, multiply_map1, multiply_map2,
	SumMap, sum_map, sum_map1, sum_map2,
	composite_jacobian,
	mul_jacobian, sum_jacobian,
	# inverse.jl
	LazyInverse, implements_inverse,
	# isomorhpism
	Isomorphism, VectorToNumber, NumberToVector,
	VectorToComplex, ComplexToVector,
	VectorToTuple, TupleToVector,
	NestedToFlat, FlatToNested,
	# jacobian.jl
	LazyJacobian, DeterminantMap, determinantmap, jacdet,
	AbsMap, absmap,
	LazyDiffVolume, diffvolume,
	NumberLike,
	to_matrix, to_vector, zeromatrix, zerovector, identitymatrix,
	# lazy.jl
	LazyMap, CompositeLazyMap, SimpleLazyMap,
	supermap,
	DerivedMap, WrappedMap,
	# map.jl
	Map, TypedMap, EuclideanMap, VectorMap,
	domaintype, codomaintype, prectype, numtype,
	convert_domaintype, convert_codomaintype, convert_prectype, convert_numtype,
	promote_map_point_pair,
	promote_and_apply, promote_maps,
	isvectorvalued_type, isvectorvalued,
	issquaremap, isoverdetermined, isunderdetermined,
	isrealmap,
	is_scalar_to_vector, is_scalar_to_scalar,
	is_vector_to_scalar, is_vector_to_vector,
	isequalmap1, isequalmap2, default_isequalmap,
	map_stencil, map_stencil_broadcast,
	# product.jl
	ProductMap, VcatMapElement,
	compatibleproductdims,
	tointernalpoint, toexternalpoint,
	productmap, productmap1, productmap2,
	VcatMap, size_as_matrix, toexternalmatrix,
	VectorProductMap, TupleProductMap,
	# affine.jl
	AbstractAffineMap,
	unsafe_matrix, unsafe_vector,
	affinematrix, affinevector,
	applymap!,
	islinearmap, isaffinemap,
	LinearMap, GenericLinearMap, ScalarLinearMap,
	VectorLinearMap, StaticLinearMap,
	Translation, ScalarTranslation, StaticTranslation,
	VectorTranslation, GenericTranslation,
	AffineMap, GenericAffineMap, ScalarAffineMap,
	VectorAffineMap, StaticAffineMap,
	# arithmetics.jl
	affine_composition,
	interval_map, bounded_interval_map,
	# basic.jl
	IdentityMap, isidentitymap,
	StaticIdentityMap, DynamicIdentityMap,
	EuclideanIdentityMap, VectorIdentityMap,
	ConstantMap, isconstantmap, mapconstant,
	ZeroMap, UnityMap, FixedConstantMap,
	# coordinates.jl
	CartToPolarMap, PolarToCartMap,
	UnitCircleMap, AngleMap, UnitDiskMap,
	# common.jl
	isrealtype,
	StaticTypes,
	hashrec,
	euclideandimension,
	convert_eltype, promotable_eltypes,
	prectype, convert_prectype, to_prectype, promote_prectype,
	numtype, convert_numtype, to_numtype, promote_numtype,
	convert_fromcartesian, convert_tocartesian,
	matrix_pinv,
	factors, nfactors, factor

import FunctionMaps:
	convert_eltype,
	convert_prectype,
	convert_numtype,
	prectype,
	numtype,
	factors,
	tointernalpoint, toexternalpoint, compatibleproductdims

# Deprecations for symbols that were renamed in FunctionMaps.jl
const AbstractMap = Map
@deprecate isreal(::Type{T}) where {T} isrealtype(T)
@deprecate isreal(m::Map) isrealmap(m)
@deprecate islinear(m::Map) islinearmap(m)
@deprecate isaffine(m::Map) isaffinemap(m)
@deprecate isidentity(m::Map) isidentitymap(m)
@deprecate isconstant(m::Map) isconstantmap(m)
@deprecate constant(m::Map) mapconstant(m)
@deprecate matrix(m::Map) affinematrix(m)
@deprecate vector(m::Map) affinevector(m)
@deprecate mapdim(m::Map) mapsize(m, 2)
