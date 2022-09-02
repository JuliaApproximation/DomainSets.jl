using Documenter, DomainSets

makedocs(
	doctest = false,
	clean = true,
	format = Documenter.HTML(),
	sitename = "DomainSets.jl",
	authors = "volunteers wanted",
	pages = Any[
			"Home" => "index.md"
	]
)


deploydocs(
    repo   = "github.com/JuliaApproximation/DomainSets.jl.git"
)
