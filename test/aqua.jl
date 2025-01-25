using Aqua

# Aqua.test_all(DomainSets)

# A few ambiguities remain in deprecated function calls. Tests to be activated
# in next breaking release when those deprecations are removed:
Aqua.test_ambiguities(DomainSets)

# The unbound tests annoyingly flag several valid uses of NTuple{N,T}
# Aqua.test_unbound_args(DomainSets)

Aqua.test_undefined_exports(DomainSets)
Aqua.test_project_extras(DomainSets)
Aqua.test_stale_deps(DomainSets)
Aqua.test_deps_compat(DomainSets)
# Most remaining piracies have to do with Domain being defined elsewhere
# Aqua.test_piracies(DomainSets)
Aqua.test_persistent_tasks(DomainSets)
