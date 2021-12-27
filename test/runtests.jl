using Test
using SafeTestsets

include("generate_example_tests.jl")

# Note: comment outer @testset to stop after first @safetestset failure
@time @testset verbose=true "Krotov.jl Package" begin

    print("\n* Example 1 (examples/simple_state_to_state.jl):")
    @time @safetestset "Example 1" begin include(joinpath("examples", "simple_state_to_state.jl")) end

    print("\n* Example 2 (examples/state_to_state_rwa.jl):")
    @time @safetestset "Example 2" begin include(joinpath("examples", "state_to_state_rwa.jl")) end

    print("\n* Example 3 (examples/rho_3states.jl):")
    if Sys.isapple()
        println("\nSkipped (macOS)")
    else
        @time @safetestset "Example 3" begin include(joinpath("examples", "rho_3states.jl")) end
    end

    print("\n* Example 4 (examples/state_to_state_parametrizations.jl):")
    @time @safetestset "Example 4" begin include(joinpath("examples", "state_to_state_parametrizations.jl")) end

    println("")

end
