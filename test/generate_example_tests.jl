import Literate
using Base.Filesystem: cp, basename

println("Start generating tests for Literate.jl examples")

EXAMPLEDIR = joinpath(@__DIR__, "..", "examples")
GENERATEDDIR = joinpath(@__DIR__, "examples")
mkpath(GENERATEDDIR)
for name in readdir(EXAMPLEDIR)
    path = abspath(joinpath(EXAMPLEDIR, name))
    if endswith(path, ".jl")
        example = path
        script = Literate.script(example, GENERATEDDIR)
    elseif isdir(path)
        folder = path
        target = joinpath(GENERATEDDIR, basename(folder))
        cp(folder, target, force=true)
        @info "cp $folder -> $target"
    end
end

println("Finished generating tests for Literate.jl examples")
