# generate examples
import Literate
using Base.Filesystem: cp, basename

println("Start generating Literate.jl examples")

EXAMPLEDIR = joinpath(@__DIR__, "..", "examples")
GENERATEDDIR = joinpath(@__DIR__, "src", "examples")
mkpath(GENERATEDDIR)
for name in readdir(EXAMPLEDIR)
    path = abspath(joinpath(EXAMPLEDIR, name))
    if endswith(path, ".jl")
        example = path
        script = Literate.script(example, GENERATEDDIR)
        code = strip(read(script, String))
        mdpost(str) = replace(str, "@__CODE__" => code)
        Literate.markdown(example, GENERATEDDIR, postprocess=mdpost)
        Literate.notebook(example, GENERATEDDIR, execute=true)
    elseif any(endswith.(path, [".png", ".jpg", ".gif"]))
        cp(path, joinpath(GENERATEDDIR, name); force=true)
    elseif isdir(path)
        folder = path
        target = joinpath(GENERATEDDIR, basename(folder))
        cp(folder, target, force=true)
    else
        @warn "ignoring $name"
    end
end

# remove any .vtu files in the generated dir (should not be deployed)
cd(GENERATEDDIR) do
    foreach(file -> endswith(file, ".vtu") && rm(file), readdir())
    foreach(file -> endswith(file, ".pvd") && rm(file), readdir())
end

println("Finished generating Literate.jl examples")
