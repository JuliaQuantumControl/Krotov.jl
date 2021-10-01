# # Example 3: Optimization of a Dissipative Quantum Gate

#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`rho_3states.ipynb`](@__NBVIEWER_ROOT_URL__/examples/rho_3states.ipynb).
#md #
#md #     Compare this example against the [same example using the `krotov`
#md #     Python package](https://qucontrol.github.io/krotov/v1.2.1/notebooks/06_example_3states.html).

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


# This example illustrates the optimization for a quantum gate in an open
# quantum system, where the dynamics is governed by the Liouville-von Neumann
# equation.

using QuantumControl
using LinearAlgebra
using Serialization
using SparseArrays

#jl using Test


# ## The two-transmon system

# We will use internal units GHz and ns. Values in GHz contain an implicit
# factor $2 \pi$, and MHz and μs are converted to GHz and ns, respectively:

const GHz = 2π;
const MHz = 0.001GHz;
const ns = 1.0;
const μs = 1000ns;
const 𝕚 = 1im;

# This implicit factor $2 \pi$ is because frequencies ($\nu$) convert to
# energies as $E = h \nu$, but our propagation routines assume a unit $\hbar =
# 1$ for energies. Thus, the factor $h / \hbar = 2 \pi$.


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

Ωre = t -> 35MHz * QuantumControl.shapes.flattop(t; T=T, t_rise=20ns);
Ωim = t -> 0.0;

L = transmon_liouvillian(Ωre, Ωim);

tlist = collect(range(0, 400ns, length=2000));

# The guess pulse looks as follows:

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

#!jl plot_control(Ωre, tlist)


# ## Optimization objectives

# Our target gate is $\Op{O} = \sqrt{\text{iSWAP}}$:

SQRTISWAP = [1  0    0   0;
             0 1/√2 𝕚/√2 0;
             0 𝕚/√2 1/√2 0;
             0  0    0   1];

# The key idea explored in the paper is that a set of three density matrices is sufficient to track the optimization
#
# $$
# \begin{align}
# \Op{\rho}_1
#     &= \sum_{i=1}^{d} \frac{2 (d-i+1)}{d (d+1)} \ketbra{i}{i} \\
# \Op{\rho}_2
#     &= \sum_{i,j=1}^{d} \frac{1}{d} \ketbra{i}{j} \\
# \Op{\rho}_3
#     &= \sum_{i=1}^{d} \frac{1}{d} \ketbra{i}{i}
# \end{align}
# $$

# In our case, $d=4$ for a two qubit-gate, and the $\ket{i}$, $\ket{j}$ are the canonical basis states $\ket{00}$, $\ket{01}$, $\ket{10}$, $\ket{11}$

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


# The three density matrices play different roles in the optimization, and, as
# shown in the paper, convergence may improve significantly by weighing the
# states relatively to each other. For this example, we place a strong emphasis
# on the optimization $\Op{\rho}_1 \rightarrow \Op{O}^\dagger \Op{\rho}_1
# \Op{O}$, by a factor of 20. This reflects that the hardest part of the
# optimization is identifying the basis in which the gate is diagonal. We will
# be using the real-part functional ($J_{T,\text{re}}$) to evaluate the success
# of $\Op{\rho}_i \rightarrow \Op{O}\Op{\rho}_i\Op{O}^\dagger$. Because
# $\Op{\rho}_1$ and $\Op{\rho}_3$ are mixed states, the Hilbert-Schmidt overlap
# will take values smaller than one in the optimal case. To compensate, we
# divide the weights by the purity of the respective states.
#
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


# ## Dynamics under the Guess Pulse

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
    objectives[1], tlist; initial_state=reshape(ρ̂₀₀, :), method=:newton,
    observables=(pop00, pop01, pop10, pop11), storage=true
);


# ## Optimization


const problem = ControlProblem(
    objectives=objectives,
    prop_method=:newton,
    use_threads=true,
    pulse_options=IdDict(
        Ωre  => Dict(
            :lambda_a => 1.0,
            :update_shape => t -> QuantumControl.shapes.flattop(t, T=T, t_rise=20ns, func=:blackman),
        ),
        Ωim  => Dict(
            :lambda_a => 1.0,
            :update_shape => t -> QuantumControl.shapes.flattop(t, T=T, t_rise=20ns, func=:blackman),
        ),
    ),
    tlist=tlist,
    iter_stop=3,
    chi=QuantumControl.functionals.chi_re!,
    J_T=QuantumControl.functionals.J_T_re,
    check_convergence= res -> begin (
            (res.J_T < 1e-3)
            && (res.converged = true)
            && (res.message="J_T < 10⁻³")
        ) end
);


opt_result = optimize(problem, method=:krotov);
#-
opt_result

# ## Continuation of Optimization

dumpdir = joinpath(@__DIR__, "dump"); mkpath(dumpdir)
dumpfile = joinpath(dumpdir, "rho_3states_opt_result.jls")
if isfile(dumpfile)
    opt_result = deserialize(dumpfile)
else
    opt_result = optimize(problem, method=:krotov, continue_from=opt_result,
                          iter_stop=3000)
    serialize(dumpfile, opt_result)
end

# ## Optimization result
