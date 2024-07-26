using Test
using SafeTestsets
using Plots

unicodeplots()

# Note: comment outer @testset to stop after first @safetestset failure
@time @testset verbose = true "Krotov.jl Package" begin

    println("\n* TLS Optimization (test_tls_optimization.jl)")
    @time @safetestset "TLS Optimization" begin
        include("test_tls_optimization.jl")
    end

    println("\n* Pulse Optimization (test_pulse_optimization.jl)")
    @time @safetestset "Pulse Optimization" begin
        include("test_pulse_optimization.jl")
    end

    println("\n* Empty Optimization (test_empty_optimization.jl)")
    @time @safetestset "Empty Optimization" begin
        include("test_empty_optimization.jl")
    end

    println("\n* Iterations (test_iterations.jl)")
    @time @safetestset "Iterations" begin
        include("test_iterations.jl")
    end

end
nothing
