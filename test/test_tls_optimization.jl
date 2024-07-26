using Test
using QuantumControl
using QuantumPropagators: ExpProp
using QuantumControl.Functionals: J_T_sm
using GRAPE
import Krotov
using LinearAlgebra
using Printf
import IOCapture
using StaticArrays: @SMatrix, @SVector

ϵ(t) = 0.2 * QuantumControl.Shapes.flattop(t, T=5, t_rise=0.3, func=:blackman);


"""Two-level-system Hamiltonian."""
function tls_hamiltonian(Ω=1.0, ϵ=ϵ)
    σ̂_z = ComplexF64[
        1  0
        0 -1
    ]
    σ̂_x = ComplexF64[
        0  1
        1  0
    ]
    Ĥ₀ = -0.5 * Ω * σ̂_z
    Ĥ₁ = σ̂_x
    return hamiltonian(Ĥ₀, (Ĥ₁, ϵ))
end;


"""Two-level-system Hamiltonian, using StaticArrays."""
function tls_hamiltonian_static(Ω=1.0, ϵ=ϵ)
    σ̂_z = @SMatrix ComplexF64[
        1  0
        0 -1
    ]
    σ̂_x = @SMatrix ComplexF64[
        0  1
        1  0
    ]
    Ĥ₀ = -0.5 * Ω * σ̂_z
    Ĥ₁ = σ̂_x
    return hamiltonian(Ĥ₀, (Ĥ₁, ϵ))
end;


@testset "TLS" begin

    println("\n==================== TLS ===========================\n")
    H = tls_hamiltonian()
    tlist = collect(range(0, 5, length=501))
    Ψ₀ = ComplexF64[1, 0]
    Ψtgt = ComplexF64[0, 1]
    problem = ControlProblem(
        [Trajectory(Ψ₀, H, target_state=Ψtgt)],
        tlist;
        iter_stop=5,
        prop_method=ExpProp,
        J_T=J_T_sm,
        check_convergence=res -> begin
            ((res.J_T < 1e-10) && (res.converged = true) && (res.message = "J_T < 10⁻¹⁰"))
        end,
    )
    res = optimize(problem; method=Krotov)
    display(res)
    @test res.J_T < 1e-3
    @test 1.0 < maximum(abs.(res.optimized_controls[1])) < 1.2
    println("===================================================\n")

end


@testset "TLS (static)" begin

    println("\n================ TLS (static) ======================\n")
    H = tls_hamiltonian_static()
    tlist = collect(range(0, 5, length=501))
    Ψ₀ = @SVector ComplexF64[1, 0]
    Ψtgt = @SVector ComplexF64[0, 1]
    problem = ControlProblem(
        [Trajectory(Ψ₀, H, target_state=Ψtgt)],
        tlist;
        iter_stop=5,
        prop_method=ExpProp,
        J_T=J_T_sm,
        check_convergence=res -> begin
            ((res.J_T < 1e-10) && (res.converged = true) && (res.message = "J_T < 10⁻¹⁰"))
        end,
    )
    res = optimize(problem; method=Krotov)
    display(res)
    @test res.J_T < 1e-3
    @test 1.0 < maximum(abs.(res.optimized_controls[1])) < 1.2
    println("===================================================\n")

end



@testset "TLS (continue from GRAPE)" begin

    println("\n============ TLS (GRAPE continuation) ============\n")
    H = tls_hamiltonian()
    tlist = collect(range(0, 5, length=501))
    Ψ₀ = ComplexF64[1, 0]
    Ψtgt = ComplexF64[0, 1]
    problem = ControlProblem(
        [Trajectory(Ψ₀, H, target_state=Ψtgt)],
        tlist;
        iter_stop=5,
        prop_method=ExpProp,
        J_T=J_T_sm,
        check_convergence=res -> begin
            ((res.J_T < 1e-10) && (res.converged = true) && (res.message = "J_T < 10⁻¹⁰"))
        end,
    )
    res_grape = optimize(problem; method=GRAPE, iter_stop=2)
    res =
        optimize(problem; method=Krotov, continue_from=res_grape, store_iter_info=["J_T"],)
    display(res)
    @test res.J_T < 1e-5
    @test abs(res.records[1][1] - res_grape.J_T) < 1e-14
    @test length(res.records) == 4
    println("===================================================\n")

end


@testset "TLS (continue with GRAPE)" begin

    println("\n=========== TLS (continue with GRAPE) ============\n")
    H = tls_hamiltonian()
    tlist = collect(range(0, 5, length=501))
    Ψ₀ = ComplexF64[1, 0]
    Ψtgt = ComplexF64[0, 1]
    problem = ControlProblem(
        [Trajectory(Ψ₀, H, target_state=Ψtgt)],
        tlist;
        iter_stop=5,
        prop_method=ExpProp,
        J_T=J_T_sm,
        check_convergence=res -> begin
            ((res.J_T < 1e-10) && (res.converged = true) && (res.message = "J_T < 10⁻¹⁰"))
        end,
    )
    res_krotov = optimize(problem; method=Krotov, iter_stop=2)
    res =
        optimize(problem; method=GRAPE, continue_from=res_krotov, store_iter_info=["J_T"],)
    display(res)
    @test res.J_T < 1e-3
    @test length(res.records) == 4
    @test abs(res.records[1][1] - res_krotov.J_T) < 1e-14
    println("===================================================\n")

end
