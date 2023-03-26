using Test
using SafeTestsets
using Plots

unicodeplots()

include("generate_example_tests.jl")

include("download_dumps.jl")

# Note: comment outer @testset to stop after first @safetestset failure
@time @testset verbose = true "Krotov.jl Package" begin

    print("\n* Pulse Optimization (test_pulse_optimization.jl)")
    @time @safetestset "Pulse Optimization" begin
        include("test_pulse_optimization.jl")
    end

    print("\n* Emptry Optimization (test_empty_optimization.jl)")
    @time @safetestset "Empty Optimization" begin
        include("test_empty_optimization.jl")
    end

    print("\n* Example 1 (examples/simple_state_to_state.jl):")
    @time @safetestset "Example 1" begin
        include(joinpath("examples", "simple_state_to_state.jl"))
    end

    print("\n* Example 2 (examples/rho_3states.jl):")
    @time @safetestset "Example 2" begin
        include(joinpath("examples", "rho_3states.jl"))
    end

    print("\n* Example 3 (examples/state_to_state_parametrizations.jl):")
    @time @safetestset "Example 3" begin
        include(joinpath("examples", "state_to_state_parametrizations.jl"))
    end

    print("\n* Example 4 (examples/perfect_entanglers.jl):")
    @time @safetestset "Example 4" begin
        include(joinpath("examples", "perfect_entanglers.jl"))
    end

    println("")

end
