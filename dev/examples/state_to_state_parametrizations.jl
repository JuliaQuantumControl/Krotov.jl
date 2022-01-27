using DrWatson
@quickactivate "KrotovTests"

using QuantumControl
using QuantumControl.Shapes: flattop
using Krotov: SquareParametrization, TanhParametrization, TanhSqParametrization, LogisticParametrization, LogisticSqParametrization
using LinearAlgebra

using Test

using PyPlot
matplotlib.use("Agg")

function plot_symmetric_parametrization_comparison()
    fig, axs = matplotlib.pyplot.subplots(figsize=(16, 3), ncols=3)

    u_vals = collect(range(-3, 3, length=101))
    ϵ_vals = collect(range(-1, 1, length=101))

    ϵ_min = -1.0
    ϵ_max = 1.0

    axs[1].plot(u_vals, u_vals, "--", color="black")
    axs[1].plot(u_vals, TanhParametrization(ϵ_min, ϵ_max).epsilon_of_u.(u_vals), label="Tanh")
    axs[1].plot(u_vals, LogisticParametrization(ϵ_min, ϵ_max).epsilon_of_u.(u_vals), label="Logistic(k=1)")
    axs[1].plot(u_vals, LogisticParametrization(ϵ_min, ϵ_max, k=4).epsilon_of_u.(u_vals), label="Logistic(k=4)")
    axs[1].set_ylim(-1.2, 1.2)
    axs[1].set_xlabel("u")
    axs[1].set_ylabel("ϵ")

    axs[2].plot(ϵ_vals, ϵ_vals, "--", color="black")
    axs[2].plot(ϵ_vals, TanhParametrization(ϵ_min, ϵ_max).u_of_epsilon.(ϵ_vals), label="Tanh")
    axs[2].plot(ϵ_vals, LogisticParametrization(ϵ_min, ϵ_max).u_of_epsilon.(ϵ_vals), label="Logistic(k=1)")
    axs[2].plot(ϵ_vals, LogisticParametrization(ϵ_min, ϵ_max, k=4).u_of_epsilon.(ϵ_vals), label="Logistic(k=4)")
    axs[2].set_ylim(-3, 3)
    axs[2].set_ylabel("u")
    axs[2].set_xlabel("ϵ")
    axs[2].legend()

    axs[3].plot(u_vals, [1.0 for _ in u_vals], "--", color="black")
    axs[3].plot(u_vals, TanhParametrization(ϵ_min, ϵ_max).de_du_derivative.(u_vals), label="Tanh")
    axs[3].plot(u_vals, LogisticParametrization(ϵ_min, ϵ_max).de_du_derivative.(u_vals), label="Logistic(k=1)")
    axs[3].plot(u_vals, LogisticParametrization(ϵ_min, ϵ_max, k=4).de_du_derivative.(u_vals), label="Logistic(k=4)")
    axs[3].set_ylim(0, 2)
    axs[3].set_xlabel("u")
    axs[3].set_ylabel("∂ϵ/∂u")

    return fig
end

function plot_positive_parametrization_comparison()
    fig, axs = matplotlib.pyplot.subplots(figsize=(16, 3), ncols=3)

    u_vals = collect(range(-3, 3, length=101))
    ϵ_vals = collect(range(0, 1, length=101))
    ϵ_max = 1.0

    axs[1].plot(u_vals, abs.(u_vals), "--", color="black")
    axs[1].plot(u_vals, TanhSqParametrization(ϵ_max).epsilon_of_u.(u_vals), label="TanhSq")
    axs[1].plot(u_vals, LogisticSqParametrization(ϵ_max).epsilon_of_u.(u_vals), label="LogisticSq(k=1)")
    axs[1].plot(u_vals, LogisticSqParametrization(ϵ_max, k=4.0).epsilon_of_u.(u_vals), label="LogisticSq(k=4)")
    axs[1].plot(u_vals, SquareParametrization().epsilon_of_u.(u_vals), label="Square")
    axs[1].set_ylim(0, 1.2)
    axs[1].set_xlabel("u")
    axs[1].set_ylabel("ϵ")

    axs[2].plot(ϵ_vals, ϵ_vals, "--", color="black")
    axs[2].plot(ϵ_vals, TanhSqParametrization(ϵ_max).u_of_epsilon.(ϵ_vals), label="TanhSq")
    axs[2].plot(ϵ_vals, LogisticSqParametrization(ϵ_max).u_of_epsilon.(ϵ_vals), label="LogisticSq(k=1)")
    axs[2].plot(ϵ_vals, LogisticSqParametrization(ϵ_max, k=4.0).u_of_epsilon.(ϵ_vals), label="LogisticSq(k=4)")
    axs[2].plot(ϵ_vals, SquareParametrization().u_of_epsilon.(ϵ_vals), label="Square")
    axs[2].set_ylim(0, 3)
    axs[2].set_ylabel("u")
    axs[2].set_xlabel("ϵ")
    axs[2].legend()

    axs[3].plot(u_vals, sign.(u_vals), "--", color="black")
    axs[3].plot(u_vals, TanhSqParametrization(ϵ_max).de_du_derivative.(u_vals), label="TanhSq")
    axs[3].plot(u_vals, LogisticSqParametrization(ϵ_max).de_du_derivative.(u_vals), label="LogisticSq(k=1)")
    axs[3].plot(u_vals, LogisticSqParametrization(ϵ_max, k=4.0).de_du_derivative.(u_vals), label="LogisticSq(k=4)")
    axs[3].plot(u_vals, SquareParametrization().de_du_derivative.(u_vals), label="Square")
    axs[3].set_ylim(-2, 2)
    axs[3].set_xlabel("u")
    axs[3].set_ylabel("∂ϵ/∂u")

    return fig
end

ϵ(t) = 0.2 * flattop(t, T=5, t_rise=0.3, func=:blackman);

"""Two-level-system Hamiltonian."""
function hamiltonian(Ω=1.0, ϵ=ϵ)
    σ̂_z = ComplexF64[1 0; 0 -1];
    σ̂_x = ComplexF64[0 1; 1  0];
    Ĥ₀ = -0.5 * Ω * σ̂_z
    Ĥ₁ = σ̂_x
    return (Ĥ₀, (Ĥ₁, ϵ))
end
;

H = hamiltonian();
@test length(H) == 2

tlist = collect(range(0, 5, length=500));

function plot_control(pulse::Vector, tlist)
    fig, ax = matplotlib.pyplot.subplots(figsize=(6, 3))
    ax.plot(tlist, pulse)
    ax.set_xlabel("time")
    ax.set_ylabel("amplitude")
    return fig
end

plot_control(ϵ::T, tlist) where T<:Function =
    plot_control([ϵ(t) for t in tlist], tlist)

function ket(label)
    result = Dict(
        "0" => Vector{ComplexF64}([1, 0]),
        "1" => Vector{ComplexF64}([0, 1]),
    )
    return result[string(label)]
end
;

@test dot(ket(0), ket(1)) ≈ 0

objectives = [
    Objective(initial_state=ket(0), generator=H, target_state=ket(1))
]

@test length(objectives) == 1

problem = ControlProblem(
    objectives=objectives,
    pulse_options=IdDict(
        ϵ  => Dict(
            :lambda_a => 5,
            :update_shape => t -> flattop(t, T=5, t_rise=0.3, func=:blackman),
            :parametrization => SquareParametrization(),
        )
    ),
    tlist=tlist,
    iter_stop=50,
    chi=QuantumControl.Functionals.chi_ss!,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence= res -> begin (
            (res.J_T < 1e-3)
            && (res.converged = true)
            && (res.message="J_T < 10⁻³")
        ) end
);

opt_result_positive = optimize(problem, method=:krotov);

opt_result_positive

@test minimum(opt_result_positive.optimized_controls[1]) ≥ 0.0
@test minimum(opt_result_positive.optimized_controls[1]) < 1e-16
@test maximum(opt_result_positive.optimized_controls[1]) > 0.0

problem_tanhsq = ControlProblem(
    objectives=objectives,
    pulse_options=IdDict(
        ϵ  => Dict(
            :lambda_a => 10,
            :update_shape => t -> flattop(t, T=5, t_rise=0.3, func=:blackman),
            :parametrization => TanhSqParametrization(3),
        )
    ),
    tlist=tlist,
    iter_stop=50,
    chi=QuantumControl.Functionals.chi_ss!,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence= res -> begin (
            (res.J_T < 1e-3)
            && (res.converged = true)
            && (res.message="J_T < 10⁻³")
        ) end
);

opt_result_tanhsq = optimize(problem_tanhsq, method=:krotov);

opt_result_tanhsq

@test minimum(opt_result_tanhsq.optimized_controls[1]) ≥ 0.0
@test minimum(opt_result_tanhsq.optimized_controls[1]) < 1e-16
@test maximum(opt_result_tanhsq.optimized_controls[1]) > 0.0
@test maximum(opt_result_tanhsq.optimized_controls[1]) < 3.0

problem_logisticsq = ControlProblem(
    objectives=objectives,
    pulse_options=IdDict(
        ϵ  => Dict(
            :lambda_a => 1,
            :update_shape => t -> flattop(t, T=5, t_rise=0.3, func=:blackman),
            :parametrization => LogisticSqParametrization(3, k=1.0),
        )
    ),
    tlist=tlist,
    iter_stop=50,
    chi=QuantumControl.Functionals.chi_ss!,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence= res -> begin (
            (res.J_T < 1e-3)
            && (res.converged = true)
            && (res.message="J_T < 10⁻³")
        ) end
);

opt_result_logisticsq = optimize(problem_logisticsq, method=:krotov);

@test minimum(opt_result_logisticsq.optimized_controls[1]) ≥ 0.0
@test minimum(opt_result_logisticsq.optimized_controls[1]) < 1e-16
@test maximum(opt_result_logisticsq.optimized_controls[1]) > 0.0
@test maximum(opt_result_logisticsq.optimized_controls[1]) < 3.0

problem_tanh = ControlProblem(
    objectives=objectives,
    pulse_options=IdDict(
        ϵ  => Dict(
            :lambda_a => 1,
            :update_shape => t -> flattop(t, T=5, t_rise=0.3, func=:blackman),
            :parametrization => TanhParametrization(-0.5, 0.5),
        )
    ),
    tlist=tlist,
    iter_stop=50,
    chi=QuantumControl.Functionals.chi_ss!,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence= res -> begin (
            (res.J_T < 1e-3)
            && (res.converged = true)
            && (res.message="J_T < 10⁻³")
        ) end
);

opt_result_tanh = optimize(problem_tanh, method=:krotov);

@test minimum(opt_result_tanh.optimized_controls[1]) > -0.5
@test maximum(opt_result_tanh.optimized_controls[1]) < 0.5

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

