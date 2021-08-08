using Krotov
using Documenter

# Generate examples
include("generate.jl")

DocMeta.setdocmeta!(Krotov, :DocTestSetup, :(using Krotov); recursive=true)

makedocs(;
    modules=[Krotov],
    authors="Michael Goerz <goerz@stanford.edu> and contributors",
    repo="https://github.com/goerz/Krotov.jl/blob/{commit}{path}#{line}",
    sitename="Krotov.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://goerz.github.io/Krotov.jl",
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

deploydocs(;
    repo="github.com/goerz/Krotov.jl",
)
