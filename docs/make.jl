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
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://QuantumControl-jl.github.io/Krotov.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Overview" => "overview.md",
        hide("Examples" => "examples/index.md",
             ["examples/simple_state_to_state.md",
              "examples/state_to_state_rwa.md"]),
        "API" => "api.md",
    ],
)

println("Finished makedocs")

rm(joinpath(@__DIR__, "build", "examples", ".gitignore"))

deploydocs(;
    repo="github.com/QuantumControl-jl/Krotov.jl",
)
