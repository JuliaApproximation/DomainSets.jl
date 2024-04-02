using Documenter
using DomainSets
using DomainSets.FunctionMaps

DocMeta.setdocmeta!(DomainSets, :DocTestSetup, :(using DomainSets); recursive=true)

makedocs(;
    modules=[DomainSets,DomainSets.FunctionMaps],
    authors="Daan Huybrechs <daan.huybrechs@kuleuven.be> and contributors",
    sitename="DomainSets.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Euclidean geometry" => "geometry.md",
        "Set operations" => "setoperations.md",
        "Public API Reference" => "api.md",
        "Internal API Reference" => "internal.md",
        "FunctionMaps.jl" => "maps.md"
    ],
)

deploydocs(
    repo   = "github.com/JuliaApproximation/DomainSets.jl.git",
    devbranch = "master"
)
