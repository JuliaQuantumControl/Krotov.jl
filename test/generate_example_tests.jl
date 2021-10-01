import Literate
import LibGit2

println("Start generating tests for Literate.jl examples")

EXAMPLEDIR = joinpath(@__DIR__, "..", "examples")
GENERATEDDIR = joinpath(@__DIR__, "examples")
mkpath(GENERATEDDIR)
for example in readdir(EXAMPLEDIR)
    if endswith(example, ".jl")
        input = abspath(joinpath(EXAMPLEDIR, example))
        script = Literate.script(input, GENERATEDDIR)
    end
end
dumpdata_dir = joinpath(GENERATEDDIR, "dump")
example_dumpdata_dir = joinpath(EXAMPLEDIR, "dump")
if !isdir(dumpdata_dir)
    if isdir(example_dumpdata_dir)
        @info "Copy $example_dumpdata_dir -> $dumpdata_dir"
        cp(example_dumpdata_dir, dumpdata_dir; force=true)
    else
        repo_url = "https://github.com/JuliaQuantumControl/Krotov.jl.git"
        @info "Clone $repo_url#data-dump -> $dumpdata_dir"
        # Direct cloning of the dump data should enable "upstream" testing in
        # CI, # i.e. `Pkg.test("Krotov")` as part of the test suite of another
        # package, without having to worry about cloning the data there
        LibGit2.clone(repo_url, dumpdata_dir, branch="data-dump")
    end
end

println("Finished generating tests for Literate.jl examples")
