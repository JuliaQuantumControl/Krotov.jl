using Test
using SafeTestsets
using Plots

unicodeplots()

include("generate_example_tests.jl")

include("download_dumps.jl")

# Note: comment outer @testset to stop after first @safetestset failure
@time @testset verbose = true "Krotov.jl Package" begin

    println("\n* Pulse Optimization (test_pulse_optimization.jl)")
    @time @safetestset "Pulse Optimization" begin
        include("test_pulse_optimization.jl")
    end

    println("\n* Empty Optimization (test_empty_optimization.jl)")
    @time @safetestset "Empty Optimization" begin
        include("test_empty_optimization.jl")
    end

    println("\n* Example 1 (examples/simple_state_to_state.jl):")
    @time @safetestset "Example 1 (simple_state_to_state)" begin
        include(joinpath("examples", "simple_state_to_state.jl"))
    end

    println("\n* Example 2 (examples/rho_3states.jl):")
    @time @safetestset "Example 2 (rho_3states)" begin
        include(joinpath("examples", "rho_3states.jl"))
    end

    println("\n* Example 3 (examples/state_to_state_parametrizations.jl):")
    @time @safetestset "Example 3 (state_to_state_parametrization)" begin
        include(joinpath("examples", "state_to_state_parametrizations.jl"))
    end

    println("\n* Example 4 (examples/perfect_entanglers.jl):")
    @time @safetestset "Example 4 (perfect_entanglers)" begin
        include(joinpath("examples", "perfect_entanglers.jl"))
    end

    println("")

end
