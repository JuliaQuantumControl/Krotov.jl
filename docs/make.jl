using QuantumControl
using QuantumPropagators
using Krotov
using Documenter
using Documenter.HTMLWriter: KaTeX
using DocumenterCitations
using DocumenterInterLinks
using Pkg
using Plots

gr()
ENV["GKSwstype"] = "100"


PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
VERSION = PROJECT_TOML["version"]
NAME = PROJECT_TOML["name"]
AUTHORS = join(PROJECT_TOML["authors"], ", ") * " and contributors"
GITHUB = "https://github.com/JuliaQuantumControl/Krotov.jl"

DEV_OR_STABLE = "stable/"
if endswith(VERSION, "dev")
    DEV_OR_STABLE = "dev/"
end


function org_inv(pkgname)
    objects_inv =
        joinpath(@__DIR__, "..", "..", "$pkgname.jl", "docs", "build", "objects.inv")
    if isfile(objects_inv)
        return ("https://juliaquantumcontrol.github.io/$pkgname.jl/dev/", objects_inv,)
    else
        return "https://juliaquantumcontrol.github.io/$pkgname.jl/$DEV_OR_STABLE"
    end
end

links = InterLinks(
    "Julia" => (
        "https://docs.julialang.org/en/v1/",
        "https://docs.julialang.org/en/v1/objects.inv",
        joinpath(@__DIR__, "src", "inventories", "Julia.toml"),
    ),
    "TimerOutputs" => (
        "https://github.com/KristofferC/TimerOutputs.jl",
        joinpath(@__DIR__, "src", "inventories", "TimerOutputs.toml")
    ),
    "QuantumPropagators" => org_inv("QuantumPropagators"),
    "QuantumControl" => org_inv("QuantumControl"),
    "GRAPE" => org_inv("GRAPE"),
    "Examples" => "https://juliaquantumcontrol.github.io/QuantumControlExamples.jl/$DEV_OR_STABLE",
    "ComponentArrays" => (
        "https://jonniedie.github.io/ComponentArrays.jl/stable/",
        "https://jonniedie.github.io/ComponentArrays.jl/stable/objects.inv",
        joinpath(@__DIR__, "src", "inventories", "ComponentArrays.toml")
    ),
    "RecursiveArrayTools" => (
        "https://docs.sciml.ai/RecursiveArrayTools/stable/",
        "https://docs.sciml.ai/RecursiveArrayTools/stable/objects.inv",
        joinpath(@__DIR__, "src", "inventories", "RecursiveArrayTools.toml")
    ),
)

println("Starting makedocs")

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)

PAGES = [
    "Home" => "index.md",
    "Overview" => "overview.md",
    "Examples" => "examples.md",
    "API" => "api.md",
    "References" => "references.md",
    hide("externals.md"),
]

makedocs(;
    plugins=[bib, links],
    modules=[Krotov],
    authors=AUTHORS,
    sitename="Krotov.jl",
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://juliaquantumcontrol.github.io/Krotov.jl",
        assets=[
            "assets/citations.css",
            asset(
                "https://juliaquantumcontrol.github.io/QuantumControl.jl/dev/assets/topbar/topbar.css"
            ),
            asset(
                "https://juliaquantumcontrol.github.io/QuantumControl.jl/dev/assets/topbar/topbar.js"
            ),
        ],
        size_threshold_ignore=["externals.md"],
        mathengine=KaTeX(
            Dict(
                :macros => Dict(
                    "\\Op" => "\\hat{#1}",
                    "\\ket" => "\\vert#1\\rangle",
                    "\\bra" => "\\langle#1\\vert",
                    "\\Im" => "\\operatorname{Im}",
                    "\\Re" => "\\operatorname{Re}",
                ),
            ),
        ),
        footer="[$NAME.jl]($GITHUB) v$VERSION docs powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).",
    ),
    pages=PAGES,
    warnonly=true,
)

println("Finished makedocs")

deploydocs(; repo="github.com/JuliaQuantumControl/Krotov.jl", devbranch="master")
