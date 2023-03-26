using QuantumControlBase
using Krotov
using Documenter
using QuantumCitations
using Pkg
using Plots

gr()
ENV["GKSwstype"] = "100"

include(joinpath("..", "test", "download_dumps.jl"))

# Generate examples
include("generate.jl")

DocMeta.setdocmeta!(Krotov, :DocTestSetup, :(using Krotov); recursive=true)

PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
VERSION = PROJECT_TOML["version"]
NAME = PROJECT_TOML["name"]
AUTHORS = join(PROJECT_TOML["authors"], ", ") * " and contributors"
GITHUB = "https://github.com/JuliaQuantumControl/Krotov.jl"

println("Starting makedocs")

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)

makedocs(
    bib;
    authors=AUTHORS,
    sitename="Krotov.jl",
    modules=[Krotov],
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://juliaquantumcontrol.github.io/Krotov.jl",
        assets=String["assets/citations.css"],
        footer="[$NAME.jl]($GITHUB) v$VERSION docs powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).",
        mathengine=KaTeX()
    ),
    pages=[
        "Home" => "index.md",
        "Overview" => "overview.md",
        "Examples" => [
            "List of Examples" => "examples/index.md",
            "Example 1 (TLS)" => "examples/simple_state_to_state.md",
            "Example 2 (Diss. Gate)" => "examples/rho_3states.md",
            "Example 3 (Parametrization)" => "examples/state_to_state_parametrizations.md",
            "Example 4 (PE)" => "examples/perfect_entanglers.md",
        ],
        "API" => "api.md",
        "References" => "references.md",
    ]
)

println("Finished makedocs")

rm(joinpath(@__DIR__, "build", "examples", ".gitignore"))

deploydocs(; repo="github.com/JuliaQuantumControl/Krotov.jl")
