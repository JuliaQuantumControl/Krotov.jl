using QuantumPropagators
using QuantumControlBase
using Krotov
using LinearAlgebra
using SparseArrays

using Test

const GHz = 2π;
const MHz = 0.001GHz;
const ns = 1.0;
const μs = 1000ns;
const 𝕚 = 1im;

function transmon_liouvillian(Ωre, Ωim;
        N=5,             # number of qubit levels
        ω₁=4.3796GHz,    # qubit frequency 1
        ω₂=4.6137GHz,    # qubit frequency 2
        ωd=4.4985GHz,    # drive frequency
        δ₁=-239.3MHz,    # anharmonicity 1
        δ₂=-242.8MHz,    # anharmonicity 2
        J=-2.3MHz,       # effective qubit-qubit coupling
        γ₁₁=(1/38.0μs),  # decay rate for qubit 1
        γ₁₂=(1/32.0μs),  # decay rate for qubit 2
        γ₂₁=(1/29.5μs),  # dephasing rate for qubit 1
        γ₂₂=(1/16.0μs),  # dephasing time for qubit 2
    )

    ⊗(A, B) = kron(A, B)
    𝟙 = SparseMatrixCSC{ComplexF64, Int64}(sparse(I, N, N))

    b̂₁ = spdiagm(1 => complex.(sqrt.(collect(1:N-1)))) ⊗ 𝟙
    b̂₂ = 𝟙 ⊗ spdiagm(1 => complex.(sqrt.(collect(1:N-1))))
    b̂₁⁺ = sparse(b̂₁'); b̂₂⁺ = sparse(b̂₂')
    n̂₁ = sparse(b̂₁' * b̂₁); n̂₂ = sparse(b̂₂' * b̂₂)
    n̂₁² = sparse(n̂₁ * n̂₁); n̂₂² = sparse(n̂₂ * n̂₂)
    b̂₁⁺_b̂₂ = sparse(b̂₁' * b̂₂); b̂₁_b̂₂⁺ = sparse(b̂₁ * b̂₂')

    Ĥ₀ = sparse(
        (ω₁ - ωd - δ₁/2) * n̂₁ + (δ₁/2) * n̂₁²
        + (ω₂ - ωd - δ₂/2) * n̂₂ + (δ₂/2) * n̂₂²
        + J * (b̂₁⁺_b̂₂ + b̂₁_b̂₂⁺)
    )

    Ĥ₁re = (1/2) * (b̂₁ + b̂₁⁺ + b̂₂ + b̂₂⁺)
    Ĥ₁im = (𝕚/2) * (b̂₁⁺ - b̂₁ + b̂₂⁺ - b̂₂)

    H = (Ĥ₀, (Ĥ₁re, Ωre), (Ĥ₁im, Ωim))

    c_ops = [√γ₁₁ * b̂₁, √γ₁₂ * b̂₂, √γ₂₁ * n̂₁, √γ₂₂ * n̂₂]

    return liouvillian(H, c_ops; convention=:TDSE)

end

const T = 400ns;

Ωre = t -> 35MHz * flattop(t; T=T, t_rise=20ns);
Ωim = t -> 0.0;

L = transmon_liouvillian(Ωre, Ωim);

tlist = collect(range(0, 400ns, length=2000));

using PyPlot
matplotlib.use("Agg")

function plot_control(pulse::Vector, tlist)
    fig, ax = matplotlib.pyplot.subplots(figsize=(6, 3))
    ax.plot(tlist, pulse)
    ax.set_xlabel("time")
    ax.set_ylabel("amplitude")
    return fig
end

plot_control(ϵ::T, tlist) where T<:Function =
    plot_control([ϵ(t) for t in tlist], tlist)

SQRTISWAP = [1  0    0   0;
             0 1/√2 𝕚/√2 0;
             0 𝕚/√2 1/√2 0;
             0  0    0   1];

function ket(i::Int64; N=5)
    Ψ = zeros(ComplexF64, N)
    Ψ[i+1] = 1
    return Ψ
end;

ket(i::Int64, j::Int64; N=5) = kron(ket(i; N=N), ket(j; N=N));

bra(args...; N=5) = adjoint(ket(args..., N=N));

const basis_labels = [(0, 0), (0, 1), (1, 0), (1, 1)];
const basis = [ket(labels...) for labels in basis_labels];
const d = length(basis);

const basis_tgt = [sum([SQRTISWAP[i,j] * basis[i] for i ∈ 1:d]) for j ∈ 1:d];


const ρ̂₁ = sum([(2*(d-i+1)/(d*(d+1))) * basis[i] * adjoint(basis[i]) for i ∈ 1:d]);
const ρ̂₂ = sum([(1/d) * basis[i] * adjoint(basis[j]) for i ∈ 1:d for j ∈ 1:d]);
const ρ̂₃ = sum([(1/d) * basis[i] * adjoint(basis[i]) for i ∈ 1:d]);

const ρ̂₁_tgt = sum([(2*(d-i+1)/(d*(d+1))) * basis_tgt[i] * adjoint(basis_tgt[i]) for i ∈ 1:d]);
const ρ̂₂_tgt = sum([(1/d) * basis_tgt[i] * adjoint(basis_tgt[j]) for i ∈ 1:d for j ∈ 1:d]);
const ρ̂₃_tgt = sum([(1/d) * basis_tgt[i] * adjoint(basis_tgt[i]) for i ∈ 1:d]);

weights = Float64[20, 1, 1];
weights *= length(weights) / sum(weights); # manual normalization
weights ./= [0.3, 1.0, 0.25]; # purities

const objectives = [
    WeightedObjective(
        initial_state=reshape(ρ̂₁,:),
        generator=L,
        target_state=reshape(ρ̂₁_tgt,:),
        weight=weights[1]
    ),
    WeightedObjective(
        initial_state=reshape(ρ̂₂,:),
        generator=L,
        target_state=reshape(ρ̂₂_tgt,:),
        weight=weights[2]
    ),
    WeightedObjective(
        initial_state=reshape(ρ̂₃,:),
        generator=L,
        target_state=reshape(ρ̂₃_tgt,:),
        weight=weights[3]
    )
];

ρ̂₀₀ = ket(0, 0) * adjoint(ket(0, 0));
ρ̂₀₁ = ket(0, 1) * adjoint(ket(0, 1));
ρ̂₁₀ = ket(1, 0) * adjoint(ket(1, 0));
ρ̂₁₁ = ket(1, 1) * adjoint(ket(1, 1));

function as_matrix(ρ⃗)
    N = isqrt(length(ρ⃗))
    return reshape(ρ⃗, N, N)
end;

pop00 = ρ⃗ -> real(tr(as_matrix(ρ⃗) * ρ̂₀₀));
pop01 = ρ⃗ -> real(tr(as_matrix(ρ⃗) * ρ̂₀₁));
pop10 = ρ⃗ -> real(tr(as_matrix(ρ⃗) * ρ̂₁₀));
pop11 = ρ⃗ -> real(tr(as_matrix(ρ⃗) * ρ̂₁₁));


rho_00_expvals = propagate(
    reshape(ρ̂₀₀, :), obj_genfunc(objectives[1], tlist), tlist; method=:newton,
    observables=(pop00, pop01, pop10, pop11), storage=true
);

const problem = ControlProblem(
    objectives=objectives,
    prop_method=:newton,
    use_threads=true,
    pulse_options=IdDict(
        Ωre  => Dict(
            :lambda_a => 1.0,
            :update_shape => t -> flattop(t, T=T, t_rise=20ns, func=:blackman),
        ),
        Ωim  => Dict(
            :lambda_a => 1.0,
            :update_shape => t -> flattop(t, T=T, t_rise=20ns, func=:blackman),
        ),
    ),
    tlist=tlist,
    iter_stop=3,
    chi=chi_re!,
    J_T=J_T_re,
    check_convergence= res -> begin (
            (res.J_T < 1e-3)
            && (res.converged = true)
            && (res.message="J_T < 10⁻³")
        ) end
);


opt_result = optimize_pulses(problem);

opt_result

using Serialization
dumpdir = joinpath(@__DIR__, "dump"); mkpath(dumpdir)
serialize(joinpath(dumpdir, "./opt_result.jls"), opt_result)
opt_result_prev = deserialize(joinpath(dumpdir, "./opt_result.jls"))

opt_result2 = optimize_pulses(problem, continue_from=opt_result_prev, iter_stop=5);

opt_result2

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

