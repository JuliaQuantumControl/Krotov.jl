using Krotov
using Test

@testset verbose=true "Krotov.jl" begin

    @testset verbose=true "Examples" begin
        include("test_examples.jl")
    end

end
