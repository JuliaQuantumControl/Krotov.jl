using Test
using SafeTestsets

# Note: comment outer @testset to stop after first @safetestset failure
@time @testset verbose=true "Krotov.jl Package" begin

    # Note: this test is copied from QuantumPropagators.jl. It works with the
    # master version of QuantumPropagators.jl, but not with v0.0.1. This test
    # is only here in order to detect whether we're successfully testing
    # against the master version of QuantumPropagators.jl. It can be removed
    # once that has been worked out.
    print("\n* Propagation (test_prop.jl):")
    @time @safetestset "Propagattion" begin include("test_prop.jl") end

    print("\n* Examples (test_examples.jl):")
    @time @safetestset "Examples" begin include("test_examples.jl") end

end
