using Test
using SafeTestsets

include("generate_example_tests.jl")

# Note: comment outer @testset to stop after first @safetestset failure
@time @testset verbose=true "Krotov.jl Package" begin

    # Note: this test is copied from QuantumPropagators.jl. It works with the
    # master version of QuantumPropagators.jl, but not with v0.0.1. This test
    # is only here in order to detect whether we're successfully testing
    # against the master version of QuantumPropagators.jl. It can be removed
    # once that has been worked out.
    print("\n* Propagation (test_prop.jl):")
    @time @safetestset "Propagation" begin include("test_prop.jl") end

    print("\n* Example 1 (examples/simple_state_to_state.jl):")
    @time @safetestset "Example 1" begin include(joinpath("examples", "simple_state_to_state.jl")) end

    print("\n* Example 2 (examples/state_to_state_rwa.jl):")
    @time @safetestset "Example 2" begin include(joinpath("examples", "state_to_state_rwa.jl")) end

    println("")

end
