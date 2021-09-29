import Literate

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
cp(joinpath(EXAMPLEDIR, "dump"), joinpath(GENERATEDDIR, "dump"); force=true)

println("Finished generating tests for Literate.jl examples")
