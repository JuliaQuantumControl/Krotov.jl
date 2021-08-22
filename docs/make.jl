using Krotov
using Documenter

# Generate examples
include("generate.jl")

DocMeta.setdocmeta!(Krotov, :DocTestSetup, :(using Krotov); recursive=true)

println("Starting makedocs")

makedocs(;
    modules=[Krotov],
    authors="Michael Goerz <mail@michaelgoerz.net> and contributors",
    repo="https://github.com/QuantumControl-jl/Krotov.jl/blob/{commit}{path}#{line}",
    sitename="Krotov.jl",
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://QuantumControl-jl.github.io/Krotov.jl",
        assets=String[],
        mathengine=KaTeX(),
    ),
    pages=[
        "Home" => "index.md",
        "Overview" => "overview.md",
        "Examples" => [
            "List of Examples" => "examples/index.md",
            "Example 1 (TLS)" => "examples/simple_state_to_state.md",
            "Example 2 (RWA)" => "examples/state_to_state_rwa.md"],
        "API" => "api.md",
    ],
)

println("Finished makedocs")

rm(joinpath(@__DIR__, "build", "examples", ".gitignore"))

deploydocs(;
    repo="github.com/QuantumControl-jl/Krotov.jl",
)
