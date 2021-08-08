# Test the scripts
import Literate

EXAMPLEDIR = joinpath(@__DIR__, "..", "examples")
GENERATEDDIR = joinpath(@__DIR__, "examples")
mkpath(GENERATEDDIR)
for example in readdir(EXAMPLEDIR)
    if endswith(example, ".jl")
        input = abspath(joinpath(EXAMPLEDIR, example))
        script = Literate.script(input, GENERATEDDIR)
    end
end


module TestSimpleStateToState
    using Test
    mktempdir() do dir
        cd(dir) do
            include(joinpath(@__DIR__, "examples", "simple_state_to_state.jl"))
        end
    end
end

module TestStateToStateRWA
    using Test
    mktempdir() do dir
        cd(dir) do
            include(joinpath(@__DIR__, "examples", "state_to_state_rwa.jl"))
        end
    end
end
