using DrWatson
@quickactivate "KrotovTests"
using Test; println("")

const GHz = 2π
const MHz = 0.001GHz
const ns = 1.0
const μs = 1000ns;

⊗ = kron
const 𝕚 = 1im
const N = 6  # levels per transmon

using LinearAlgebra
using SparseArrays
using QuantumControl


function transmon_hamiltonian(;
    Ωre,
    Ωim,
    N=N,  # levels per transmon
    ω₁=4.380GHz,
    ω₂=4.614GHz,
    ωd=4.498GHz,
    α₁=-210MHz,
    α₂=-215MHz,
    J=-3MHz,
    λ=1.03,
    use_sparse=:auto
)
    𝟙 = SparseMatrixCSC{ComplexF64,Int64}(sparse(I, N, N))
    b̂₁ = spdiagm(1 => complex.(sqrt.(collect(1:N-1)))) ⊗ 𝟙
    b̂₂ = 𝟙 ⊗ spdiagm(1 => complex.(sqrt.(collect(1:N-1))))
    b̂₁⁺ = sparse(b̂₁')
    b̂₂⁺ = sparse(b̂₂')
    n̂₁ = sparse(b̂₁' * b̂₁)
    n̂₂ = sparse(b̂₂' * b̂₂)
    n̂₁² = sparse(n̂₁ * n̂₁)
    n̂₂² = sparse(n̂₂ * n̂₂)
    b̂₁⁺_b̂₂ = sparse(b̂₁' * b̂₂)
    b̂₁_b̂₂⁺ = sparse(b̂₁ * b̂₂')

    ω̃₁ = ω₁ - ωd
    ω̃₂ = ω₂ - ωd

    Ĥ₀ = sparse(
        (ω̃₁ - α₁ / 2) * n̂₁ +
        (α₁ / 2) * n̂₁² +
        (ω̃₂ - α₂ / 2) * n̂₂ +
        (α₂ / 2) * n̂₂² +
        J * (b̂₁⁺_b̂₂ + b̂₁_b̂₂⁺)
    )

    Ĥ₁re = (1 / 2) * (b̂₁ + b̂₁⁺ + λ * b̂₂ + λ * b̂₂⁺)
    Ĥ₁im = (𝕚 / 2) * (b̂₁⁺ - b̂₁ + λ * b̂₂⁺ - λ * b̂₂)

    if ((N < 5) && (use_sparse ≢ true)) || use_sparse ≡ false
        H = hamiltonian(Array(Ĥ₀), (Array(Ĥ₁re), Ωre), (Array(Ĥ₁im), Ωim))
    else
        H = hamiltonian(Ĥ₀, (Ĥ₁re, Ωre), (Ĥ₁im, Ωim))
    end
    return H

end;

using QuantumControl.Shapes: flattop

function guess_pulses(; T=400ns, E₀=35MHz, dt=0.1ns, t_rise=15ns)

    tlist = collect(range(0, T, step=dt))
    Ωre = t -> E₀ * flattop(t, T=T, t_rise=t_rise)
    Ωim = t -> 0.0

    return tlist, Ωre, Ωim

end

tlist, Ωre_guess, Ωim_guess = guess_pulses();

using Plots
Plots.default(
    linewidth               = 3,
    size                    = (550, 300),
    legend                  = :right,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.8),
)

function plot_complex_pulse(tlist, Ω; time_unit=:ns, ampl_unit=:MHz, kwargs...)

    ax1 = plot(
        tlist ./ eval(time_unit),
        abs.(Ω) ./ eval(ampl_unit);
        label="|Ω|",
        xlabel="time ($time_unit)",
        ylabel="amplitude ($ampl_unit)",
        kwargs...
    )

    ax2 = plot(
        tlist ./ eval(time_unit),
        angle.(Ω) ./ π;
        label="ϕ(Ω)",
        xlabel="time ($time_unit)",
        ylabel="phase (π)"
    )

    plot(ax1, ax2, layout=(2, 1))

end

plot_complex_pulse(tlist, Ωre_guess.(tlist) + 𝕚 * Ωim_guess.(tlist))

function ket(i::Int64; N=N)
    Ψ = zeros(ComplexF64, N)
    Ψ[i+1] = 1
    return Ψ
end

function ket(indices::Int64...; N=N)
    Ψ = ket(indices[1]; N=N)
    for i in indices[2:end]
        Ψ = Ψ ⊗ ket(i; N=N)
    end
    return Ψ
end

function ket(label::AbstractString; N=N)
    indices = [parse(Int64, digit) for digit in label]
    return ket(indices...; N=N)
end

basis = [ket("00"), ket("01"), ket("10"), ket("11")];

H = transmon_hamiltonian(Ωre=Ωre_guess, Ωim=Ωim_guess);

objectives = [Objective(; initial_state=Ψ, generator=H) for Ψ ∈ basis];

using QuantumControl.WeylChamber: D_PE, gate_concurrence, unitarity
using QuantumControl.Functionals: gate_functional

J_T_PE = gate_functional(D_PE; unitarity_weight=0.5);

using QuantumControl: propagate_objectives

guess_states = propagate_objectives(objectives, tlist; use_threads=true);

U_guess = [basis[i] ⋅ guess_states[j] for i = 1:4, j = 1:4]

gate_concurrence(U_guess)

@test gate_concurrence(U_guess) < 0.9

using QuantumControl.WeylChamber: weyl_chamber_region
weyl_chamber_region(U_guess)

@test weyl_chamber_region(U_guess) == "W0*"

1 - unitarity(U_guess)

@test round(1 - unitarity(U_guess), digits=1) ≈ 0.1

D_PE(U_guess)

0.5 * D_PE(U_guess) + 0.5 * (1 - unitarity(U_guess))

J_T_PE(guess_states, objectives)

@test 0.4 < J_T_PE(guess_states, objectives) < 0.5
@test 0.5 * D_PE(U_guess) + 0.5 * (1-unitarity(U_guess)) ≈ J_T_PE(guess_states, objectives) atol=1e-15

problem = ControlProblem(
    objectives=objectives,
    lambda_a=10.0,
    update_shape=(t -> flattop(t, T=400ns, t_rise=15ns, func=:blackman)),
    tlist=tlist,
    iter_stop=100,
    J_T=J_T_PE,
    check_convergence=res -> begin
        (
            (res.J_T > res.J_T_prev) &&
            (res.converged = true) &&
            (res.message = "Loss of monotonic convergence")
        )
        (
            (res.J_T <= 1e-3) &&
            (res.converged = true) &&
            (res.message = "Found a perfect entangler")
        )
    end,
    use_threads=true,
);

opt_result = @optimize_or_load(datadir("PE_OCT.jld2"), problem; method=:Krotov);

opt_result

Ω_opt = opt_result.optimized_controls[1] + 𝕚 * opt_result.optimized_controls[2]

plot_complex_pulse(tlist, Ω_opt)

using QuantumControl.Controls: substitute, get_controls

opt_states = propagate_objectives(
    substitute(
        objectives,
        IdDict(zip(get_controls(objectives), opt_result.optimized_controls))
    ),
    tlist;
    use_threads=true
);

U_opt = [basis[i] ⋅ opt_states[j] for i = 1:4, j = 1:4];

gate_concurrence(U_opt)
@test round(gate_concurrence(U_opt), digits=3) ≈ 1.0

1 - unitarity(U_opt)

@test round(1 - unitarity(U_opt), digits=2) ≈ 0.01

J_T_C = U -> 0.5 * (1 - gate_concurrence(U)) + 0.5 * (1 - unitarity(U));

using QuantumControl.Functionals: make_gate_chi

chi_C = make_gate_chi(J_T_C, objectives);

opt_result_direct = @optimize_or_load(
    datadir("PE_OCT_direct.jld2"),
    problem;
    method=:Krotov,
    J_T=gate_functional(J_T_C),
    chi=chi_C
);

opt_result_direct

opt_states_direct = propagate_objectives(
    substitute(
        objectives,
        IdDict(zip(get_controls(objectives), opt_result_direct.optimized_controls))
    ),
    tlist;
    use_threads=true
);

U_opt_direct = [basis[i] ⋅ opt_states_direct[j] for i = 1:4, j = 1:4];

gate_concurrence(U_opt_direct)
@test round(gate_concurrence(U_opt_direct), digits=3) ≈ 1.0

1 - unitarity(U_opt_direct)
@test round(1 - unitarity(U_opt_direct), digits=3) ≈ 0.002

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

