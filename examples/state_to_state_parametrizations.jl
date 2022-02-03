# # Example 3: Pulse Parametrization

#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`simple_state_to_state.ipynb`](@__NBVIEWER_ROOT_URL__/examples/state_to_state_parametrization.ipynb).

#md # ``\gdef\op#1{\hat{#1}}``
#md # ``\gdef\init{\text{init}}``
#md # ``\gdef\tgt{\text{tgt}}``

#nb # $
#nb # \newcommand{tr}[0]{\operatorname{tr}}
#nb # \newcommand{diag}[0]{\operatorname{diag}}
#nb # \newcommand{abs}[0]{\operatorname{abs}}
#nb # \newcommand{pop}[0]{\operatorname{pop}}
#nb # \newcommand{aux}[0]{\text{aux}}
#nb # \newcommand{opt}[0]{\text{opt}}
#nb # \newcommand{tgt}[0]{\text{tgt}}
#nb # \newcommand{init}[0]{\text{init}}
#nb # \newcommand{lab}[0]{\text{lab}}
#nb # \newcommand{rwa}[0]{\text{rwa}}
#nb # \newcommand{bra}[1]{\langle#1\vert}
#nb # \newcommand{ket}[1]{\vert#1\rangle}
#nb # \newcommand{Bra}[1]{\left\langle#1\right\vert}
#nb # \newcommand{Ket}[1]{\left\vert#1\right\rangle}
#nb # \newcommand{Braket}[2]{\left\langle #1\vphantom{#2}\mid{#2}\vphantom{#1}\right\rangle}
#nb # \newcommand{op}[1]{\hat{#1}}
#nb # \newcommand{Op}[1]{\hat{#1}}
#nb # \newcommand{dd}[0]{\,\text{d}}
#nb # \newcommand{Liouville}[0]{\mathcal{L}}
#nb # \newcommand{DynMap}[0]{\mathcal{E}}
#nb # \newcommand{identity}[0]{\mathbf{1}}
#nb # \newcommand{Norm}[1]{\lVert#1\rVert}
#nb # \newcommand{Abs}[1]{\left\vert#1\right\vert}
#nb # \newcommand{avg}[1]{\langle#1\rangle}
#nb # \newcommand{Avg}[1]{\left\langle#1\right\rangle}
#nb # \newcommand{AbsSq}[1]{\left\vert#1\right\vert^2}
#nb # \newcommand{Re}[0]{\operatorname{Re}}
#nb # \newcommand{Im}[0]{\operatorname{Im}}
#nb # $


# This example illustrates the parametrization of control pulses as a
# form of amplitude constraint.

using DrWatson
@quickactivate "KrotovTests"
#-
using QuantumControl
using QuantumControl.Shapes: flattop
using Krotov: SquareParametrization, TanhParametrization, TanhSqParametrization, LogisticParametrization, LogisticSqParametrization
using LinearAlgebra

#jl using Test; println("")
using Plots

# ## Parametrizations

# ## Symmetric Bounded Controls


function plot_symmetric_parametrization_comparison()

    u_vals = collect(range(-3, 3, length=101))
    ϵ_vals = collect(range(-1, 1, length=101))

    ϵ_min = -1.0
    ϵ_max = 1.0

    legend_args = Dict(
        :legend => :topleft,
        :foreground_color_legend => nothing,
        :background_color_legend => RGBA(1, 1, 1, 0.8)
    )

    pnl1 = plot(u_vals, u_vals; linestyle=:dash, color="black", label="", xlabel="u", ylabel="ϵ", legend=false)
    plot!(pnl1, u_vals, TanhParametrization(ϵ_min, ϵ_max).epsilon_of_u.(u_vals), label="Tanh")
    plot!(pnl1,u_vals, LogisticParametrization(ϵ_min, ϵ_max).epsilon_of_u.(u_vals), label="Logistic(k=1)")
    plot!(pnl1, u_vals, LogisticParametrization(ϵ_min, ϵ_max, k=4).epsilon_of_u.(u_vals), label="Logistic(k=4)")
    ylims!(pnl1, (-1.2, 1.2))

    pnl2 = plot(ϵ_vals, ϵ_vals; linestyle=:dash, color="black", label="", xlabel="ϵ", ylabel="u", legend_args...)
    plot!(pnl2, ϵ_vals, TanhParametrization(ϵ_min, ϵ_max).u_of_epsilon.(ϵ_vals), label="Tanh")
    plot!(pnl2, ϵ_vals, LogisticParametrization(ϵ_min, ϵ_max).u_of_epsilon.(ϵ_vals), label="Logistic(k=1)")
    plot!(pnl2, ϵ_vals, LogisticParametrization(ϵ_min, ϵ_max, k=4).u_of_epsilon.(ϵ_vals), label="Logistic(k=4)")
    ylims!(pnl2, (-3, 3))

    pnl3 = plot(u_vals, [1.0 for _ in u_vals]; linestyle=:dash, color="black", label="", xlabel="u", ylabel="∂ϵ/∂u", legend=false)
    plot!(pnl3, u_vals, TanhParametrization(ϵ_min, ϵ_max).de_du_derivative.(u_vals), label="Tanh")
    plot!(pnl3, u_vals, LogisticParametrization(ϵ_min, ϵ_max).de_du_derivative.(u_vals), label="Logistic(k=1)")
    plot!(pnl3, u_vals, LogisticParametrization(ϵ_min, ϵ_max, k=4).de_du_derivative.(u_vals), label="Logistic(k=4)")
    ylims!(pnl3, (0, 2))

    plot(pnl1, pnl2, pnl3, layout=(1, 3), size = (1000, 300),
         left_margin=20Plots.px, bottom_margin=20Plots.px)

end
#-
fig = plot_symmetric_parametrization_comparison()
#jl display(fig)

# ## Positive (Bounded) Controls

function plot_positive_parametrization_comparison()

    u_vals = collect(range(-3, 3, length=101))
    ϵ_vals = collect(range(0, 1, length=101))
    ϵ_max = 1.0

    legend_args = Dict(
        :legend => :topleft,
        :foreground_color_legend => nothing,
        :background_color_legend => RGBA(1, 1, 1, 0.8)
    )

    pnl1 = plot(u_vals, abs.(u_vals); linestyle=:dash, color="black", label="", xlabel="u", ylabel="ϵ", legend=false)
    plot!(pnl1, u_vals, TanhSqParametrization(ϵ_max).epsilon_of_u.(u_vals), label="TanhSq")
    plot!(pnl1, u_vals, LogisticSqParametrization(ϵ_max).epsilon_of_u.(u_vals), label="LogisticSq(k=1)")
    plot!(pnl1, u_vals, LogisticSqParametrization(ϵ_max, k=4.0).epsilon_of_u.(u_vals), label="LogisticSq(k=4)")
    plot!(pnl1, u_vals, SquareParametrization().epsilon_of_u.(u_vals), label="Square")
    ylims!(pnl1, (0, 1.2))

    pnl2 = plot(ϵ_vals, ϵ_vals; linestyle=:dash, color="black", label="", xlabel="ϵ", ylabel="u", legend_args...)
    plot!(pnl2, ϵ_vals, TanhSqParametrization(ϵ_max).u_of_epsilon.(ϵ_vals), label="TanhSq")
    plot!(pnl2, ϵ_vals, LogisticSqParametrization(ϵ_max).u_of_epsilon.(ϵ_vals), label="LogisticSq(k=1)")
    plot!(pnl2, ϵ_vals, LogisticSqParametrization(ϵ_max, k=4.0).u_of_epsilon.(ϵ_vals), label="LogisticSq(k=4)")
    plot!(pnl2, ϵ_vals, SquareParametrization().u_of_epsilon.(ϵ_vals), label="Square")
    ylims!(pnl2, (0, 3))

    pnl3 = plot(u_vals, sign.(u_vals); linestyle=:dash, color="black", label="", xlabel="u", ylabel="∂ϵ/∂u", legend=false)
    plot!(pnl3, u_vals, TanhSqParametrization(ϵ_max).de_du_derivative.(u_vals), label="TanhSq")
    plot!(pnl3, u_vals, LogisticSqParametrization(ϵ_max).de_du_derivative.(u_vals), label="LogisticSq(k=1)")
    plot!(pnl3, u_vals, LogisticSqParametrization(ϵ_max, k=4.0).de_du_derivative.(u_vals), label="LogisticSq(k=4)")
    plot!(pnl3, u_vals, SquareParametrization().de_du_derivative.(u_vals), label="Square")
    ylims!(pnl3, (-2, 2))

    plot(pnl1, pnl2, pnl3, layout=(1, 3), size = (1000, 300),
         left_margin=20Plots.px, bottom_margin=20Plots.px)

end
#-
fig = plot_positive_parametrization_comparison()
#jl display(fig)

# ## Two-level Hamiltonian

# We consider the Hamiltonian $\op{H}_{0} = - \frac{\omega}{2} \op{\sigma}_{z}$, representing
# a simple qubit with energy level splitting $\omega$ in the basis
# $\{\ket{0},\ket{1}\}$. The control field $\epsilon(t)$ is assumed to couple via
# the Hamiltonian $\op{H}_{1}(t) = \epsilon(t) \op{\sigma}_{x}$ to the qubit,
# i.e., the control field effectively drives transitions between both qubit
# states.
#
# We we will use

ϵ(t) = 0.2 * flattop(t, T=5, t_rise=0.3, func=:blackman);


#-
"""Two-level-system Hamiltonian."""
function hamiltonian(Ω=1.0, ϵ=ϵ)
    σ̂_z = ComplexF64[1 0; 0 -1];
    σ̂_x = ComplexF64[0 1; 1  0];
    Ĥ₀ = -0.5 * Ω * σ̂_z
    Ĥ₁ = σ̂_x
    return (Ĥ₀, (Ĥ₁, ϵ))
end
;
#-

H = hamiltonian();
#jl @test length(H) == 2

# The control field here switches on from zero at $t=0$ to it's maximum amplitude
# 0.2 within the time period 0.3 (the switch-on shape is half a [Blackman pulse](https://en.wikipedia.org/wiki/Window_function#Blackman_window)).
# It switches off again in the time period 0.3 before the
# final time $T=5$). We use a time grid with 500 time steps between 0 and $T$:

tlist = collect(range(0, 5, length=500));

#-
function plot_control(pulse::Vector, tlist)
    plot(tlist, pulse, xlabel="time", ylabel="amplitude", legend=false)
end

plot_control(ϵ::T, tlist) where T<:Function =
    plot_control([ϵ(t) for t in tlist], tlist)

plot_control(H[2][2], tlist)

# ## Optimization target

# The `krotov` package requires the goal of the optimization to be described by a
# list of `Objective` instances. In this example, there is only a single
# objective: the state-to-state transfer from initial state $\ket{\Psi_{\init}} =
# \ket{0}$ to the target state $\ket{\Psi_{\tgt}} = \ket{1}$, under the dynamics
# of the Hamiltonian $\op{H}(t)$:

function ket(label)
    result = Dict(
        "0" => Vector{ComplexF64}([1, 0]),
        "1" => Vector{ComplexF64}([0, 1]),
    )
    return result[string(label)]
end
;

#-
#jl @test dot(ket(0), ket(1)) ≈ 0
#-

objectives = [
    Objective(initial_state=ket(0), generator=H, target_state=ket(1))
]

#-
#jl @test length(objectives) == 1
#-


# ## Square-parametrization for positive pulses


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
#-
opt_result_positive = optimize(problem, method=:krotov);
#-
opt_result_positive

# We can plot the optimized field:

#-
#!jl plot_control(opt_result_positive.optimized_controls[1], tlist)
#-

#-
#jl @test minimum(opt_result_positive.optimized_controls[1]) ≥ 0.0
#jl @test minimum(opt_result_positive.optimized_controls[1]) < 1e-16
#jl @test maximum(opt_result_positive.optimized_controls[1]) > 0.0
#-

# ## Tanh-Square-Parametrization for positive amplitude-constrained pulses


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
#-
opt_result_tanhsq = optimize(problem_tanhsq, method=:krotov);
#-
opt_result_tanhsq

# We can plot the optimized field:

#-
#!jl plot_control(opt_result_tanhsq.optimized_controls[1], tlist)
#-

#-
#jl @test minimum(opt_result_tanhsq.optimized_controls[1]) ≥ 0.0
#jl @test minimum(opt_result_tanhsq.optimized_controls[1]) < 1e-16
#jl @test maximum(opt_result_tanhsq.optimized_controls[1]) > 0.0
#jl @test maximum(opt_result_tanhsq.optimized_controls[1]) < 3.0
#-

# ## Logistic-Square-Parametrization for positive amplitude-constrained pulses

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
#-
opt_result_logisticsq = optimize(problem_logisticsq, method=:krotov);
# We can plot the optimized field:

#-
#!jl plot_control(opt_result_logisticsq.optimized_controls[1], tlist)
#-

#-
#jl @test minimum(opt_result_logisticsq.optimized_controls[1]) ≥ 0.0
#jl @test minimum(opt_result_logisticsq.optimized_controls[1]) < 1e-16
#jl @test maximum(opt_result_logisticsq.optimized_controls[1]) > 0.0
#jl @test maximum(opt_result_logisticsq.optimized_controls[1]) < 3.0
#-
# ## Tanh-parametrization for amplitude-constrained pulses

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
#-
opt_result_tanh = optimize(problem_tanh, method=:krotov);
#-
#!jl plot_control(opt_result_tanh.optimized_controls[1], tlist)
#-
#jl @test minimum(opt_result_tanh.optimized_controls[1]) > -0.5
#jl @test maximum(opt_result_tanh.optimized_controls[1]) < 0.5
#-
