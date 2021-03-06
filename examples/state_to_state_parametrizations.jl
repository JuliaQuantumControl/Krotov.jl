# # Example 3: Pulse Parametrization

#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`state_to_state_parametrizations.ipynb`](@__NBVIEWER_ROOT_URL__/examples/state_to_state_parametrizations.ipynb).

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
using Krotov:
    SquareParametrization,
    TanhParametrization,
    TanhSqParametrization,
    LogisticParametrization,
    LogisticSqParametrization
using LinearAlgebra

#-
#jl using Test; println("")
using Plots
Plots.default(
    linewidth               = 3,
    size                    = (550, 300),
    legend                  = :top,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.8),
)
#-

# ## Parametrizations

# ## Symmetric Bounded Controls

#-
include(joinpath(@__DIR__, "plots", "symmetric_parametrization_comparison.jl"))  # hide
fig = plot_symmetric_parametrization_comparison()  # hide
#jl display(fig)

# ## Positive (Bounded) Controls

#-
include(joinpath(@__DIR__, "plots", "positive_parametrization_comparison.jl"))  # hide
fig = plot_positive_parametrization_comparison()  # hide
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

??(t) = 0.2 * flattop(t, T=5, t_rise=0.3, func=:blackman);


#-
"""Two-level-system Hamiltonian."""
function hamiltonian(??=1.0, ??=??)
    ????_z = ComplexF64[
        1  0
        0 -1
    ]
    ????_x = ComplexF64[
        0  1
        1  0
    ]
    H????? = -0.5 * ?? * ????_z
    H????? = ????_x
    return (H?????, (H?????, ??))
end;
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

plot_control(??::T, tlist) where {T<:Function} = plot_control([??(t) for t in tlist], tlist)

plot_control(H[2][2], tlist)

# ## Optimization target

# The `krotov` package requires the goal of the optimization to be described by a
# list of `Objective` instances. In this example, there is only a single
# objective: the state-to-state transfer from initial state $\ket{\Psi_{\init}} =
# \ket{0}$ to the target state $\ket{\Psi_{\tgt}} = \ket{1}$, under the dynamics
# of the Hamiltonian $\op{H}(t)$:

function ket(label)
    result = Dict("0" => Vector{ComplexF64}([1, 0]), "1" => Vector{ComplexF64}([0, 1]),)
    return result[string(label)]
end;

#-
#jl @test dot(ket(0), ket(1)) ??? 0
#-

objectives = [Objective(initial_state=ket(0), generator=H, target_state=ket(1))]

#-
#jl @test length(objectives) == 1
#-


# ## Square-parametrization for positive pulses


problem = ControlProblem(
    objectives=objectives,
    pulse_options=IdDict(
        ?? => Dict(
            :lambda_a => 5,
            :update_shape => t -> flattop(t, T=5, t_rise=0.3, func=:blackman),
            :parametrization => SquareParametrization(),
        )
    ),
    tlist=tlist,
    iter_stop=50,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence=res -> begin
        ((res.J_T < 1e-3) && (res.converged = true) && (res.message = "J_T < 10?????"))
    end
);
#-
opt_result_positive, file = @optimize_or_load(
    datadir(),
    problem;
    method=:krotov,
    filename="parametrization#opt_result_positive.jld2"
);
#-
opt_result_positive

# We can plot the optimized field:

#-
#!jl plot_control(opt_result_positive.optimized_controls[1], tlist)
#-

#-
#jl @test minimum(opt_result_positive.optimized_controls[1]) ??? 0.0
#jl @test minimum(opt_result_positive.optimized_controls[1]) < 1e-16
#jl @test maximum(opt_result_positive.optimized_controls[1]) > 0.0
#-

# ## Tanh-Square-Parametrization for positive amplitude-constrained pulses


problem_tanhsq = ControlProblem(
    objectives=objectives,
    pulse_options=IdDict(
        ?? => Dict(
            :lambda_a => 10,
            :update_shape => t -> flattop(t, T=5, t_rise=0.3, func=:blackman),
            :parametrization => TanhSqParametrization(3),
        )
    ),
    tlist=tlist,
    iter_stop=50,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence=res -> begin
        ((res.J_T < 1e-3) && (res.converged = true) && (res.message = "J_T < 10?????"))
    end
);
#-
opt_result_tanhsq, file = @optimize_or_load(
    datadir(),
    problem_tanhsq;
    method=:krotov,
    filename="parametrization#opt_result_tanhsq.jld2"
);
#-
opt_result_tanhsq

# We can plot the optimized field:

#-
#!jl plot_control(opt_result_tanhsq.optimized_controls[1], tlist)
#-

#-
#jl @test minimum(opt_result_tanhsq.optimized_controls[1]) ??? 0.0
#jl @test minimum(opt_result_tanhsq.optimized_controls[1]) < 1e-16
#jl @test maximum(opt_result_tanhsq.optimized_controls[1]) > 0.0
#jl @test maximum(opt_result_tanhsq.optimized_controls[1]) < 3.0
#-

# ## Logistic-Square-Parametrization for positive amplitude-constrained pulses

problem_logisticsq = ControlProblem(
    objectives=objectives,
    pulse_options=IdDict(
        ?? => Dict(
            :lambda_a => 1,
            :update_shape => t -> flattop(t, T=5, t_rise=0.3, func=:blackman),
            :parametrization => LogisticSqParametrization(3, k=1.0),
        )
    ),
    tlist=tlist,
    iter_stop=50,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence=res -> begin
        ((res.J_T < 1e-3) && (res.converged = true) && (res.message = "J_T < 10?????"))
    end
);
#-
opt_result_logisticsq, file = @optimize_or_load(
    datadir(),
    problem_logisticsq;
    method=:krotov,
    filename="parametrization#opt_result_logisticsq.jld2"
);
# We can plot the optimized field:

#-
#!jl plot_control(opt_result_logisticsq.optimized_controls[1], tlist)
#-

#-
#jl @test minimum(opt_result_logisticsq.optimized_controls[1]) ??? 0.0
#jl @test minimum(opt_result_logisticsq.optimized_controls[1]) < 1e-16
#jl @test maximum(opt_result_logisticsq.optimized_controls[1]) > 0.0
#jl @test maximum(opt_result_logisticsq.optimized_controls[1]) < 3.0
#-
# ## Tanh-parametrization for amplitude-constrained pulses

problem_tanh = ControlProblem(
    objectives=objectives,
    pulse_options=IdDict(
        ?? => Dict(
            :lambda_a => 1,
            :update_shape => t -> flattop(t, T=5, t_rise=0.3, func=:blackman),
            :parametrization => TanhParametrization(-0.5, 0.5),
        )
    ),
    tlist=tlist,
    iter_stop=50,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence=res -> begin
        ((res.J_T < 1e-3) && (res.converged = true) && (res.message = "J_T < 10?????"))
    end
);
#-
opt_result_tanh, file = @optimize_or_load(
    datadir(),
    problem_tanh;
    method=:krotov,
    filename="parametrization#opt_result_tanh.jld2"
);
#-
#!jl plot_control(opt_result_tanh.optimized_controls[1], tlist)
#-
#jl @test minimum(opt_result_tanh.optimized_controls[1]) > -0.5
#jl @test maximum(opt_result_tanh.optimized_controls[1]) < 0.5
#-
