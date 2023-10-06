const PROJECTDIR = dirname(Base.active_project())
projectdir(names...) = joinpath(PROJECTDIR, names...)
datadir(names...) = projectdir("data", names...)

using QuantumControl
using QuantumControl.Shapes: flattop
using QuantumControl.Generators
using QuantumControl.Controls
using QuantumControl.PulseParametrizations:
    SquareParametrization,
    TanhParametrization,
    TanhSqParametrization,
    LogisticParametrization,
    LogisticSqParametrization,
    ParametrizedAmplitude
using LinearAlgebra

using Test; println("")
using Plots
Plots.default(
    linewidth               = 3,
    size                    = (550, 300),
    legend                  = :top,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.8),
)

include(joinpath(@__DIR__, "plots", "symmetric_parametrization_comparison.jl"))  # hide
fig = plot_symmetric_parametrization_comparison()  # hide
display(fig)

include(joinpath(@__DIR__, "plots", "positive_parametrization_comparison.jl"))  # hide
fig = plot_positive_parametrization_comparison()  # hide
display(fig)

ϵ(t) = 0.2 * flattop(t, T=5, t_rise=0.3, func=:blackman);

"""Two-level-system Hamiltonian."""
function tls_hamiltonian(; Ω=1.0, ampl=ϵ)
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
    return hamiltonian(Ĥ₀, (Ĥ₁, ampl))
end;

H = tls_hamiltonian();
@test length(H.ops) == 2

tlist = collect(range(0, 5, length=500));

function plot_amplitude(ampl, tlist)
    plot(tlist, discretize(ampl, tlist), xlabel="time", ylabel="amplitude", legend=false)
end

plot_amplitude(ϵ, tlist)

function ket(label)
    result = Dict("0" => Vector{ComplexF64}([1, 0]), "1" => Vector{ComplexF64}([0, 1]),)
    return result[string(label)]
end;

@test dot(ket(0), ket(1)) ≈ 0

objectives = [Objective(initial_state=ket(0), generator=H, target_state=ket(1))]

@test length(objectives) == 1

a = ParametrizedAmplitude(
    ϵ,
    tlist;
    parametrization=SquareParametrization(),
    parameterize=true
)

function plot_amplitude(ampl::ParametrizedAmplitude, tlist)
    plot(
        tlist,
        discretize(Array(ampl), tlist),
        xlabel="time",
        ylabel="amplitude",
        legend=false
    )
end

plot_amplitude(a, tlist)

problem = ControlProblem(
    objectives=substitute(objectives, IdDict(ϵ => a)),
    lambda_a=5,
    update_shape=(t -> flattop(t, T=5, t_rise=0.3, func=:blackman)),
    tlist=tlist,
    iter_stop=50,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence=res -> begin
        ((res.J_T < 1e-3) && (res.converged = true) && (res.message = "J_T < 10⁻³"))
    end
);

opt_result_positive = @optimize_or_load(
    datadir("parametrization#opt_result_positive.jld2"),
    problem;
    method=:krotov
);

opt_result_positive

amplitude = Array(substitute(a, IdDict(a.control => opt_result_positive.optimized_controls[1])))
@test minimum(amplitude) ≥ 0.0
@test minimum(amplitude) < 1e-16
@test maximum(amplitude) > 0.0

a = ParametrizedAmplitude(
    ϵ,
    tlist;
    parametrization=TanhSqParametrization(3),
    parameterize=true
)

problem_tanhsq = ControlProblem(
    objectives=substitute(objectives, IdDict(ϵ => a)),
    lambda_a=10,
    update_shape=(t -> flattop(t, T=5, t_rise=0.3, func=:blackman)),
    tlist=tlist,
    iter_stop=50,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence=res -> begin
        ((res.J_T < 1e-3) && (res.converged = true) && (res.message = "J_T < 10⁻³"))
    end
);

opt_result_tanhsq = @optimize_or_load(
    datadir("parametrization#opt_result_tanhsq.jld2"),
    problem_tanhsq;
    method=:krotov
);

opt_result_tanhsq

amplitude = Array(substitute(a, IdDict(a.control => opt_result_tanhsq.optimized_controls[1])))
@test minimum(amplitude) ≥ 0.0
@test minimum(amplitude) < 1e-16
@test maximum(amplitude) > 0.0
@test maximum(opt_result_tanhsq.optimized_controls[1]) < 3.0

a = ParametrizedAmplitude(
    ϵ,
    tlist;
    parametrization=LogisticSqParametrization(3, k=1.0),
    parameterize=true
)

problem_logisticsq = ControlProblem(
    objectives=substitute(objectives, IdDict(ϵ => a)),
    lambda_a=1,
    update_shape=(t -> flattop(t, T=5, t_rise=0.3, func=:blackman)),
    tlist=tlist,
    iter_stop=50,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence=res -> begin
        ((res.J_T < 1e-3) && (res.converged = true) && (res.message = "J_T < 10⁻³"))
    end
);

opt_result_logisticsq = @optimize_or_load(
    datadir("parametrization#opt_result_logisticsq.jld2"),
    problem_logisticsq;
    method=:krotov
);

amplitude = Array(substitute(a, IdDict(a.control => opt_result_logisticsq.optimized_controls[1])))
@test minimum(amplitude) ≥ 0.0
@test minimum(amplitude) < 1e-16
@test maximum(amplitude) > 0.0
@test maximum(amplitude) < 3.0

a = ParametrizedAmplitude(
    ϵ,
    tlist;
    parametrization=TanhParametrization(-0.5, 0.5),
    parameterize=true
)

problem_tanh = ControlProblem(
    objectives=substitute(objectives, IdDict(ϵ => a)),
    lambda_a=1,
    update_shape=(t -> flattop(t, T=5, t_rise=0.3, func=:blackman)),
    tlist=tlist,
    iter_stop=50,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence=res -> begin
        ((res.J_T < 1e-3) && (res.converged = true) && (res.message = "J_T < 10⁻³"))
    end
);

opt_result_tanh = @optimize_or_load(
    datadir("parametrization#opt_result_tanh.jld2"),
    problem_tanh;
    method=:krotov
);

amplitude = Array(substitute(a, IdDict(a.control => opt_result_tanh.optimized_controls[1])))
@test minimum(amplitude) > -0.5
@test maximum(amplitude) < 0.5

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
