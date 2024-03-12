using Aqua

# Aqua.test_all(FunctionMaps)

Aqua.test_ambiguities(FunctionMaps)
# The unbound tests annoyingly flag several valid uses of NTuple{N,T}
# Aqua.test_unbound_args(FunctionMaps)
Aqua.test_undefined_exports(FunctionMaps)
Aqua.test_project_extras(FunctionMaps)
Aqua.test_stale_deps(FunctionMaps)
Aqua.test_deps_compat(FunctionMaps)
# Piracies to be removed in next breaking release
# Aqua.test_piracies(FunctionMaps)
Aqua.test_persistent_tasks(FunctionMaps)
