using QuantumPropagators
using QuantumControlBase
using Krotov

using Test

const σ̂_z = ComplexF64[1 0; 0 -1];
const σ̂_x = ComplexF64[0 1; 1  0];

"""Two-level-system Hamiltonian."""
function hamiltonian(Ω=1.0, E0=0.2)
    Ĥ₀ = -0.5 * Ω * σ̂_z
    Ĥ₁ = σ̂_x
    ϵ(t) = E0 * flattop(t, T=5, t_rise=0.3, func=:blackman)
    return (Ĥ₀, (Ĥ₁, ϵ))
end
;

H = hamiltonian();
@test length(H) == 2

tlist = collect(range(0, 5, length=500));

using PyPlot

function plot_pulse(pulse::Vector, tlist)
    fig, ax = matplotlib.pyplot.subplots(figsize=(6, 3))
    ax.plot(tlist, pulse)
    ax.set_xlabel("time")
    ax.set_ylabel("amplitude")
    return fig
end

plot_pulse(ϵ::T, tlist) where T<:Function =
    plot_pulse([ϵ(t) for t in tlist], tlist)

function ket(label)
    result = Dict(
        "0" => Vector{ComplexF64}([1, 0]),
        "1" => Vector{ComplexF64}([0, 1]),
    )
    return result[string(label)]
end
;

using LinearAlgebra
@test dot(ket(0), ket(1)) ≈ 0

objectives = [Objective(initial_state=ket(0), generator=H, target=ket(1))]

@test length(objectives) == 1

problem = ControlProblem(
    objectives=objectives,
    pulse_options=Dict(
        H[2][2]  => Dict(
            :lambda_a => 5,
            :update_shape => t -> flattop(t, T=5, t_rise=0.3, func=:blackman),
        )
    ),
    tlist=tlist,
);

guess_dynamics = propagate(
        objectives[1], problem.tlist;
        storage=true, observables=(Ψ->abs.(Ψ).^2, )
)

function plot_population(pop0::Vector, pop1::Vector, tlist)
    fig, ax = matplotlib.pyplot.subplots(figsize=(6, 3))
    ax.plot(tlist, pop0, label="0")
    ax.plot(tlist, pop1, label="1")
    ax.legend()
    ax.set_xlabel("time")
    ax.set_ylabel("population")
    return fig
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

