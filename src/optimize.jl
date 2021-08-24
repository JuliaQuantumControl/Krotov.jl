using QuantumControlBase
using QuantumPropagators: initpropwrk, propstep!, init_storage, write_to_storage!, get_from_storage!
import Base.Threads.@threads

# Krotov workspace (for internal use)
struct KrotovWrk

    # a copy of the objectives
    objectives :: Vector{QuantumControlBase.Objective}

    # the adjoint objectives, containing the adjoint generators for the
    # backward propagation
    adjoint_objectives :: Vector{QuantumControlBase.Objective}

    # time grid for propagated states 0, dt, ... T
    tlist :: Vector{Float64}

    # Tuple of the original controls (probably functions)
    controls :: Tuple

    # controls discretized on intervals of tlist
    guess_pulses :: Vector{Any}

    # The optimized pulses
    pulses :: Vector{Any}

    # map of controls to options
    pulse_options :: Dict{Any, Any}

    #################################
    # scratch objects, per objective:

    state :: Vector{Any}

    G :: Vector{Any}

    vals_dict :: Vector{Any}

    fw_storage :: Vector{Any}  # forward storage array (per objective)

    fw_storage2 :: Vector{Any}  # forward storage array (per objective)

    bw_storage :: Vector{Any}  # backward storage array (per objective)

    prop_wrk :: Vector{Any}

    function KrotovWrk(problem::QuantumControlBase.ControlProblem)
        objectives = [obj for obj in problem.objectives]
        adjoint_objectives = [adjoint(obj) for obj in problem.objectives]
        controls = getcontrols(objectives)
        tlist = problem.tlist
        guess_pulses = [
            discretize_on_midpoints(control, tlist) for control in controls
        ]
        pulses = [copy(pulse) for pulse in guess_pulses]
        pulse_options = problem.pulse_options
        state = [similar(obj.initial_state) for obj in objectives]
        zero_vals = IdDict(control => zero(guess_pulses[i][1]) for (i, control) in enumerate(controls))
        G = [setcontrolvals(obj.generator, zero_vals) for obj in objectives]
        vals_dict = [copy(zero_vals) for _ in objectives]
        # TODO: forward_storage only if g_b != 0
        fw_storage = [init_storage(obj.initial_state, tlist) for obj in objectives]
        # TODO: second forward storage only if second order
        fw_storage2 = [init_storage(obj.initial_state, tlist) for obj in objectives]
        bw_storage = [init_storage(obj.initial_state, tlist) for obj in objectives]
        prop_wrk = [initpropwrk(state[i], tlist, G[i]) for i in 1:length(objectives)]
        # TODO: Result object
        new(objectives, adjoint_objectives, tlist, controls,
            guess_pulses, pulses, pulse_options, state, G, vals_dict,
            fw_storage, fw_storage2, bw_storage, prop_wrk)
    end

end


"""Use Krotov's method to optimize the given optimization problem.

```julia
result = optimize_pulses(control_problem; chi_func, kwargs...)
```

optimizes the control problem defined in `control_problem`. The optimization
functional is given implicitly via the mandatory keyword argument `chi_func`.

The `control_problem` is set up as a
[`QuantumControlBase.ControlProblem`](@ref).

# Optional Keyword Arguments

The following keyword arguments are supported (with default values):

* `sigma=nothing`: Function that calculate the second-order contribution. If not given,
   the first-order Krotov method is used.
* `iter_start=0`: the initial iteration number
* `iter_stop=5000`: the maximum iteration number

"""
function optimize_pulses(control_problem; kwargs...)
    #=chi_constructor = kwargs[:chi_constructor]=#
    sigma = get(kwargs, :sigma, nothing)
    iter_start = get(kwargs, :iter_start, 0)
    iter_stop = get(kwargs, :iter_stop, 5000)
    skip_initial_forward_propagation = get(kwargs, :skip_initial_forward_propagation, false)

    wrk = KrotovWrk(control_problem)
    # TODO: if continuing previous optimization, ammend work from existing Result

    i = iter_start

    # TODO: tic-toc time (for infohook)
    if skip_initial_forward_propagation
        @info "Skipping initial forward propagation"
    else
        @info "Initial Forward Propagation"
        @threads for (k, obj) in collect(enumerate(wrk.objectives))
            krotov_initial_fw_prop!(wrk.state[k], obj.initial_state, k, wrk)
        end
    end

    # TODO: tau_vals
    # TODO: if sigma, fw_storage0 = fw_storage
    # TODO: set up Result
    # TODO: info_hook

    return wrk # XXX (premature exit for debugging)

    converged = (iter_stop <= i)

    @info "Optimization"
    while !converged
        # TODO: tic-toc time (for infohook)
        @debug "Krotov iteration $i"
        krotov_iteration(wrk)
        # TODO: info_hook
        # TODO: check convergence
        i = i + 1
    end

    # TODO: return result

end

# The dynamical generator for the forward propagation
function _fw_gen(k, n, wrk) # TODO
    vals_dict = wrk.vals_dict[k]
    ϵ = wrk.guess_pulses
    t = wrk.tlist
    for (l, control) in enumerate(wrk.controls)
        vals_dict[control] = ϵ[l][n]
    end
    dt = t[n+1] - t[n]
    setcontrolvals!(wrk.G[k], wrk.objectives[k].generator, vals_dict)
    return wrk.G[k], dt
end

function krotov_initial_fw_prop!(ϕₖ, ϕₖⁱⁿ, k, wrk)
    Φ₀ = wrk.fw_storage[k]
    copyto!(ϕₖ,  ϕₖⁱⁿ)
    (Φ₀ ≠ nothing) && write_to_storage!(Φ₀, 1,  ϕₖⁱⁿ)
    N_T = length(wrk.tlist) - 1
    for n = 1:N_T
        G, dt = _fw_gen(k, n, wrk)
        propstep!(ϕₖ, G, dt, wrk.prop_wrk[k])
        (Φ₀ ≠ nothing) && write_to_storage!(Φ₀, n+1, ϕₖ)
    end
end

# The dynamical generator for the backward propagation
function _bw_gen(k, n, wrk) # TODO: check if this is correct
    vals_dict = wrk.vals_dict[k]
    ϵ = wrk.guess_pulses
    t = wrk.tlist
    for (l, control) in enumerate(wrk.controls)
        vals_dict[control] = ϵ[l][n]
    end
    dt = t[n+1] - t[n]
    setcontrolvals!(wrk.G[k], wrk.adjoint_objectives[k].generator, vals_dict)
    return G, -dt
end

function krotov_iteration(wrk)

    ϕ = wrk.state  # assumed to contain the results of forward propagation
    chi_constructor = ϕ -> ϕ  # TODO: get function from arguments
    χ = chi_constructor(ϕ)
    N_T = length(wrk.tlist) - 1
    N = length(wrk.objectives)
    L = length(wrk.controls)
    X = wrk.bw_storage

    # backward propagation
    @threads for k = 1:N
        write_to_storage!(X[k], N_T+1, χ[k])
        for n = N_T:-1:1
            local (G, dt) = _bw_gen(k, n, wrk)
            propstep!(χ[k], G, dt, wrk.prop_work[k])
            write_to_storage!(X[k], n, χ[k])
        end
    end

    # pulse update and forward propagation

    @threads for k = 1:N
        copyto!(ϕ[k], wrk.objectives[k].initial_state)
    end

    for n = 1:N_T
        χ = wrk.state
        for k = 1:N
            get_from_storage!(χ[k], X[k], n)
        end
        Δϵₙ = [_update() for l in 1:L]
        # TODO: second order
        # TODO: apply pulse update
        @threads for k = 1:N
            # TODO: fw-prop and storage
        end
        # TODO: update sigma
    end
end
