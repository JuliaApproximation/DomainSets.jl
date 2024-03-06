
# In this file we list all functions defined in FunctionMaps that
# were defined in the namespace of DomainSets

using .FunctionMaps:
	MapRef,
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
	AbstractMap, Map, TypedMap, EuclideanMap, VectorMap,
	domaintype, codomaintype, prectype, numtype,
	convert_domaintype, convert_codomaintype, convert_prectype, convert_numtype,
	promote_map_point_pair,
	promote_and_apply, promote_maps,
	isvectorvalued_type, isvectorvalued,
	issquaremap, isoverdetermined, isunderdetermined,
	isreal,
	is_scalar_to_vector, is_scalar_to_scalar,
	is_vector_to_scalar, is_vector_to_vector,
	isequalmap1, isequalmap2, default_isequalmap,
	map_stencil, map_stencil_broadcast,
	# product.jl
	ProductMap, VcatMapElement,
	compatibleproductdims,
	tointernalpoint, toexternalpoint,
	productmap, productmap1, productmap2,
	VcatMap, mapdim, size_as_matrix, toexternalmatrix,
	VectorProductMap, TupleProductMap,
	# affine.jl
	AbstractAffineMap,
	unsafe_matrix, unsafe_vector, matrix, vector,
	affinematrix, affinevector,
	applymap!,
	islinear, isaffine,
	islinearmap, isaffinemap,
	LinearMap, GenericLinearMap, ScalarLinearMap,
	VectorLinearMap, StaticLinearMap,
	Translation, ScalarTranslation, StaticTranslation,
	VectorTranslation, GenericTranslation,
	AffineMap, GenericAffineMap, ScalarAffineMap,
	VectorAffineMap, StaticAffineMap,
	# arithmetics.jl
	affine_composition,
	# basic.jl
	IdentityMap, isidentity, isidentitymap,
	StaticIdentityMap, DynamicIdentityMap,
	EuclideanIdentityMap, VectorIdentityMap,
	ConstantMap, isconstantmap, mapconstant,
	isconstant, constant,
	ZeroMap, UnityMap, FixedConstantMap,
	# common.jl
	StaticTypes,
	hashrec,
	euclideandimension,
	convert_eltype, promotable_eltypes,
	prectype, convert_prectype, to_prectype, promote_prectype,
	numtype, convert_numtype, to_numtype, promote_numtype,
	convert_fromcartesian, convert_tocartesian,
	matrix_pinv,
	factors, nfactors, factor

import .FunctionMaps:
	convert_eltype,
	convert_prectype,
	convert_numtype,
	prectype,
	numtype,
	factors,
	isreal,
	isequalmap,
	mapsize, applymap,
	jacobian, jacobian!, diffvolume,
	leftinverse, rightinverse, inverse,
	tointernalpoint, toexternalpoint, compatibleproductdims
