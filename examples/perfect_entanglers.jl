# # Example 4: Optimization for a perfect entangler
#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`perfect_entanglers.ipynb`](@__NBVIEWER_ROOT_URL__/examples/perfect_entanglers.ipynb).
#md #
#md #     Compare this example against the [same example using GRAPE](https://juliaquantumcontrol.github.io/GRAPE.jl/stable/examples/perfect_entanglers/).

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

using DrWatson
@quickactivate "KrotovTests"
#jl using Test; println("")

# This example illustrates the optimization towards a perfectly entangling
# two-qubit gate for a system of two transmon qubits with a shared transmission
# line. It uses both the indirect perfect entanglers functional shown in
# Goerz et. al., Phys. Rev. A 91, 062307 (2015) and a direct maximization of
# the gate concurrence and thus demonstrates the optimization for non-analytic
# functions via the calculation of gradients with automatic differentiation.

# ## Hamiltonian and guess pulses

# We will write the Hamiltonian in units of GHz (angular frequency; the factor
# 2?? is implicit) and ns:

const GHz = 2??
const MHz = 0.001GHz
const ns = 1.0
const ??s = 1000ns;

# The Hamiltonian and parameter are taken from Goerz et. al., Phys. Rev. A 91,
# 062307 (2015)., cf. Table 1 in that Reference.

??? = kron
const ???? = 1im
const N = 6  # levels per transmon

using LinearAlgebra
using SparseArrays

function hamiltonian(;
    ??re,
    ??im,
    N=N,  # levels per transmon
    ?????=4.380GHz,
    ?????=4.614GHz,
    ??d=4.498GHz,
    ?????=-210MHz,
    ?????=-215MHz,
    J=-3MHz,
    ??=1.03,
    use_sparse=:auto
)
    ???? = SparseMatrixCSC{ComplexF64,Int64}(sparse(I, N, N))
    b????? = spdiagm(1 => complex.(sqrt.(collect(1:N-1)))) ??? ????
    b????? = ???? ??? spdiagm(1 => complex.(sqrt.(collect(1:N-1))))
    b???????? = sparse(b?????')
    b???????? = sparse(b?????')
    n????? = sparse(b?????' * b?????)
    n????? = sparse(b?????' * b?????)
    n??????? = sparse(n????? * n?????)
    n??????? = sparse(n????? * n?????)
    b????????_b????? = sparse(b?????' * b?????)
    b?????_b???????? = sparse(b????? * b?????')

    ??????? = ????? - ??d
    ??????? = ????? - ??d

    H????? = sparse(
        (??????? - ????? / 2) * n????? +
        (????? / 2) * n??????? +
        (??????? - ????? / 2) * n????? +
        (????? / 2) * n??????? +
        J * (b????????_b????? + b?????_b????????)
    )

    H?????re = (1 / 2) * (b????? + b???????? + ?? * b????? + ?? * b????????)
    H?????im = (???? / 2) * (b???????? - b????? + ?? * b???????? - ?? * b?????)

    if ((N < 5) && (use_sparse ??? true)) || use_sparse ??? false
        H = (Array(H?????), (Array(H?????re), ??re), (Array(H?????im), ??im))
    else
        H = (H?????, (H?????re, ??re), (H?????im, ??im))
    end
    return H

end;

# We choose a pulse duration of 400 ns. The guess pulse amplitude is 35 MHz,
# with a 15 ns switch-on/-off time. The Hamiltonian is written in a rotating
# frame, so in general, the control field is allowed to be complex-valued. We
# separate this into two control fields, one for the real part and one for the
# imaginary part. Initially, the imaginary part is zero, corresponding to a
# field exactly at the frequency of the rotating frame.

using QuantumControl.Shapes: flattop

function guess_pulses(; T=400ns, E???=35MHz, dt=0.1ns, t_rise=15ns)

    tlist = collect(range(0, T, step=dt))
    ??re = t -> E??? * flattop(t, T=T, t_rise=t_rise)
    ??im = t -> 0.0

    return tlist, ??re, ??im

end

tlist, ??re_guess, ??im_guess = guess_pulses();

# We can visualize this:

using Plots
Plots.default(
    linewidth               = 3,
    size                    = (550, 300),
    legend                  = :right,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.8),
)

function plot_complex_pulse(tlist, ??; time_unit=:ns, ampl_unit=:MHz, kwargs...)

    ax1 = plot(
        tlist ./ eval(time_unit),
        abs.(??) ./ eval(ampl_unit);
        label="|??|",
        xlabel="time ($time_unit)",
        ylabel="amplitude ($ampl_unit)",
        kwargs...
    )

    ax2 = plot(
        tlist ./ eval(time_unit),
        angle.(??) ./ ??;
        label="??(??)",
        xlabel="time ($time_unit)",
        ylabel="phase (??)"
    )

    plot(ax1, ax2, layout=(2, 1))

end

plot_complex_pulse(tlist, ??re_guess.(tlist) + ???? * ??im_guess.(tlist))

# ## Logical basis for two-qubit gates

# For simplicity, we will be define the qubits in the *bare* basis, i.e.
# ignoring the static coupling $J$.

function ket(i::Int64; N=N)
    ?? = zeros(ComplexF64, N)
    ??[i+1] = 1
    return ??
end

function ket(indices::Int64...; N=N)
    ?? = ket(indices[1]; N=N)
    for i in indices[2:end]
        ?? = ?? ??? ket(i; N=N)
    end
    return ??
end

function ket(label::AbstractString; N=N)
    indices = [parse(Int64, digit) for digit in label]
    return ket(indices...; N=N)
end

basis = [ket("00"), ket("01"), ket("10"), ket("11")];

# ## Defining the optimization problem

# We define the optimization with one objective for each of the four basis
# states:

using QuantumControl

H = hamiltonian(??re=??re_guess, ??im=??im_guess);

objectives = [MinimalObjective(; initial_state=??, generator=H) for ?? ??? basis];

# Note the use of `MinimalObjective` which does not have a `target_state`. This
# is because we will be optimizing for an arbitrary perfect entangler, not for
# a specific quantum gate.

# The optimization is steered by the perfect entanglers distance measure
# $D_{PE}$, that is, the geometric distance of the quantum gate obtained from
# propagating the four basis states to the polyhedron of perfect entanglers in
# the Weyl chamber. Since the logical subspace defining the qubit is embedded
# in the larger Hilbert space of the transmon, there may be loss of population
# from the logical subspace. To counter this possibility in the optimization,
# we add a unitarity measure  to $D_PE$. The two terms are added with equal
# weight.

using QuantumControl.WeylChamber: D_PE, gate_concurrence, unitarity
using QuantumControl.Functionals: gate_functional

J_T_PE = gate_functional(D_PE; unitarity_weight=0.5);

# The `gate_functional` routines used above converts the function `D_PE` that
# receives the gate $U??$ as a 4??4 matrix into a functional of the correct from
# for the `QuantumControl.optimize` routine, which is a function of the
# propagated states.

# We can check that for the guess pulse, we are not implementing a perfect
# entangler

using QuantumControl: propagate_objectives

guess_states = propagate_objectives(objectives, tlist; use_threads=true);

U_guess = [basis[i] ??? guess_states[j] for i = 1:4, j = 1:4]

gate_concurrence(U_guess)

#jl @test gate_concurrence(U_guess) < 0.9

# We find that the guess pulse produces a gate in the `W0*` region of the Weyl
# chamber:

using QuantumControl.WeylChamber: weyl_chamber_region
weyl_chamber_region(U_guess)

#jl @test weyl_chamber_region(U_guess) == "W0*"

# That is, the region of the Weyl chamber containing controlled-phase gates with
# a phase $> ??$ (Weyl chamber coordinates $c??? > ??/2$, $c??? < ??/4$).

# This in fact allows use to use the perfect entangler functional without
# modification: if the guess pulse were in the "W1" region of the Weyl chamber,
# (close to SWAP), we would have to flip its sign, or we would optimize towards
# the local equivalence class of the SWAP gate instead of towards the perfect
# of perfect entanglers. In principle, we could use a modified functional that
# takes the absolute square of the `D_PE` term, by using
#
# ```
# J_T_PE = gate_functional(D_PE; unitarity_weight=0.5, absolute_square=true)
# ```
#
# This would specifically optimize for the *surface* of the perfect
# entanglers functional.

# The guess pulse loses about 10% of population from the logical subspace:

1 - unitarity(U_guess)

#jl @test round(1 - unitarity(U_guess), digits=1) ??? 0.1

# We can also evaluate the geometric distance to the polyhedron of perfect
# entanglers in the Weyl chamber:

D_PE(U_guess)

# Together with the unitarity measure, this is the initial value of the
# optimization functional:

0.5 * D_PE(U_guess) + 0.5 * (1 - unitarity(U_guess))
#-
J_T_PE(guess_states, objectives)

#jl @test 0.4 < J_T_PE(guess_states, objectives) < 0.5
#jl @test 0.5 * D_PE(U_guess) + 0.5 * (1-unitarity(U_guess)) ??? J_T_PE(guess_states, objectives) atol=1e-15

# ## Optimization

# Now, we formulate the full control problem

problem = ControlProblem(
    objectives=objectives,
    pulse_options=IdDict(
        ??re_guess => Dict(
            :lambda_a => 10.0,
            :update_shape => t -> flattop(t, T=400ns, t_rise=15ns, func=:blackman),
        ),
        ??im_guess => Dict(
            :lambda_a => 10.0,
            :update_shape => t -> flattop(t, T=400ns, t_rise=15ns, func=:blackman),
        ),
    ),
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

# Note that we have not not given a `chi` parameter to calculate the boundary
# condition $|???????? = -???J_T/???????????|$ that Krotov's method requires. In this case,
# the Krotov.jl package will use automatic differentiation to determine the
# derivative. In principle, the perfect entanglers function has an analytical
# derivative, but it is exceedingly laborious to calculate and implement it
# (see the [source code of the `weylchamber` Python
# package](https://github.com/qucontrol/weylchamber/blob/master/src/weylchamber/perfect_entanglers.py),
# which can evaluate that derivative).


opt_result, file =
    @optimize_or_load(datadir(), problem; method=:Krotov, filename="PE_OCT.jld2");
#-
opt_result


# ## Optimization result

# We extract the optimized control field from the optimization result and plot
# it

??_opt = opt_result.optimized_controls[1] + ???? * opt_result.optimized_controls[2]

plot_complex_pulse(tlist, ??_opt)

# We then propagate the optimized control field to analyze the resulting
# quantum gate:

opt_states = propagate_objectives(
    objectives,
    tlist;
    use_threads=true,
    controls_map=IdDict(
        ??re_guess => opt_result.optimized_controls[1],
        ??im_guess => opt_result.optimized_controls[2]
    )
);

U_opt = [basis[i] ??? opt_states[j] for i = 1:4, j = 1:4];

# We find that we have achieved a perfect entangler:

gate_concurrence(U_opt)
#jl @test round(gate_concurrence(U_opt), digits=3) ??? 1.0

# Moreover, we have reduced the population loss to ??? 1%

1 - unitarity(U_opt)

#jl @test round(1 - unitarity(U_opt), digits=2) ??? 0.01


# ## Direct maximization of the gate concurrence

# In the previous optimizations, we have optimized for a perfect entangler
# indirectly via a geometric function in the Weyl chamber. The entire reason
# that perfect entangler functional was formulated is because calculating the
# gate concurrence directly involves the eigenvalues of the unitary, see
# Kraus, Cirac, Phys. Rev. A 63, 062309 (2001) and
# Childs et al., PRA 68, 052311 (2003), which are inherently non-analytic.

# However, since Krotov.jl can use automatic differentiation, this is no longer
# an insurmountable obstacle!

# We can define a functional for a given gate `U` that combines the gate
# concurrence and (as above) a unitarity measure to penalize loss of population
# from the logical subspace:

J_T_C = U -> 0.5 * (1 - gate_concurrence(U)) + 0.5 * (1 - unitarity(U));

# In the optimization, we will convert this functional to one that takes the
# propagated states as arguments (via the `gate_functional` routine).
# We can do that same for the gradient: Let Zygote (the automatic
# differentiation framework we are using) determine the gradient for `J_T_C`
# with respect to `U`, but then analytically translate that into the derivative
# with respect to the states that we need to calculate the ??-states.

using QuantumControl.Functionals: make_gate_chi

chi_C = make_gate_chi(J_T_C, objectives);

# Running this, we again are able to find a perfect entangler.

opt_result_direct, file = @optimize_or_load(
    datadir(),
    problem;
    method=:Krotov,
    J_T=gate_functional(J_T_C),
    chi=chi_C,
    filename="PE_OCT_direct.jld2"
);
#-
opt_result_direct
#-
opt_states_direct = propagate_objectives(
    objectives,
    tlist;
    use_threads=true,
    controls_map=IdDict(
        ??re_guess => opt_result_direct.optimized_controls[1],
        ??im_guess => opt_result_direct.optimized_controls[2]
    )
);

U_opt_direct = [basis[i] ??? opt_states_direct[j] for i = 1:4, j = 1:4];
#-
gate_concurrence(U_opt_direct)
#jl @test round(gate_concurrence(U_opt_direct), digits=3) ??? 1.0
#-
1 - unitarity(U_opt_direct)
#jl @test round(1 - unitarity(U_opt_direct), digits=3) ??? 0.002
