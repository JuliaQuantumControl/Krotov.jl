using Test
using LinearAlgebra
using QuantumPropagators

@testset "random Hermitian state" begin

    SERIALIZATION = 0
    SERIALIZE_NUMPY = true
    # A normal test run is with SERIALIZATION = 0. This uses a random
    # Hamiltonian and initial state
    # Set NORMALIZATION = 1 to write the random Hamiltonian and initial state
    # to disk, and then NORMALIZATION = 2 to load it from disk (ensuring a
    # reproducible test run.
    # In combination with NORMALIZATION = 1, if SERIALIZE_NUMPY is true, the
    # Hamiltonian and initial state will also be dumped into a numpy (.npy)
    # file, which then can be read into the test_newton.py script.

    precision = 1e-10

    # Input
    N = 1000

    dt = 0.5

    if SERIALIZATION < 2
        X = rand(ComplexF64, (N, N))
        H = Hermitian(X)
        Ψ₀ = rand(ComplexF64, N)
        Ψ₀ ./= norm(Ψ₀)
    end

    @test norm(Ψ₀) ≈ 1

    # Expected
    U = exp(-im * H * dt)
    @test norm(U * U'  - one(U)) < precision  # U is unitary
    Ψ_out_expected = U * Ψ₀
    @test norm(Ψ_out_expected) ≈ 1

    Ψ = copy(Ψ₀)
    wrk = NewtonWrk(Ψ₀, 5)
    newton!(Ψ, H, dt, wrk)
    Ψ_out = copy(Ψ)
    @test norm(Ψ_out) ≈ 1

    # Comparison
    @test norm(Ψ_out - Ψ_out_expected) < precision

end
