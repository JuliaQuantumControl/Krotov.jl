using QuantumPropagators: propstep!, init_storage, write_to_storage!, get_from_storage!
using QuantumControlBase
using QuantumControlBase.ConditionalThreads: @threadsif
using LinearAlgebra
using Dates
using Printf

"""Result object returned by [`optimize_pulses`](@ref)."""
mutable struct KrotovResult{STST}
    tlist :: Vector{Float64}
    iter_start :: Int64  # the starting iteration number
    iter_stop :: Int64 # the maximum iteration number
    iter :: Int64  # the current iteration number
    secs :: Float64  # seconds that the last iteration took
    tau_vals :: Vector{ComplexF64}
    J_T :: Float64  # the current value of the final-time functional J_T
    J_T_prev :: Float64  # previous value of J_T
    guess_controls :: Vector{Vector{Float64}}
    optimized_controls :: Vector{Vector{Float64}}
    all_pulses :: Vector{Any} # TODO
    states :: Vector{STST}
    start_local_time :: DateTime
    end_local_time :: DateTime
    records :: Vector{Tuple}  # storage for info_hook to write data into at each iteration
    converged :: Bool
    message :: String

    function KrotovResult(problem)
        tlist = problem.tlist
        controls = getcontrols(problem.objectives)
        iter_start = get(problem.kwargs, :iter_start, 0)
        iter_stop = get(problem.kwargs, :iter_stop, 5000)
        iter = 0
        secs = 0
        tau_vals = Vector{ComplexF64}()
        guess_controls = [
            discretize(control, tlist) for control in controls
        ]
        J_T = 0.0
        J_T_prev = 0.0
        optimized_controls = [copy(guess) for guess in guess_controls]
        all_pulses = Vector{Any}()
        states = [similar(obj.initial_state) for obj in problem.objectives]
        start_local_time = now()
        end_local_time = now()
        records = Vector{Tuple}()
        converged = false
        message = "in progress"
        new{eltype(states)}(
            tlist, iter_start, iter_stop, iter, secs, tau_vals, J_T, J_T_prev,
            guess_controls, optimized_controls, all_pulses, states,
            start_local_time, end_local_time, records, converged, message)
    end
end

Base.show(io::IO, r::KrotovResult) = print(io, "KrotovResult<$(r.message)>")
Base.show(io::IO, ::MIME"text/plain", r::KrotovResult) = print(io, """
Krotov Optimization Result
--------------------------
- Started at $(r.start_local_time)
- Number of objectives: $(length(r.states))
- Number of iterations: $(max(r.iter - r.iter_start, 0))
- Reason for termination: $(r.message)
- Ended at $(r.end_local_time) ($(r.end_local_time - r.start_local_time))""")


# Krotov workspace (for internal use)
struct KrotovWrk{
        OT<:QuantumControlBase.AbstractControlObjective,
        AOT<:QuantumControlBase.AbstractControlObjective,
        KWT, CTRST<:Tuple, POT<:AbstractDict, STST, VDT, STORT, PRWT, GT
    }

    # a copy of the objectives
    objectives :: Vector{OT}

    # the adjoint objectives, containing the adjoint generators for the
    # backward propagation
    adjoint_objectives :: Vector{AOT}

    # The kwargs from the control problem
    kwargs :: KWT

    # Tuple of the original controls (probably functions)
    controls :: CTRST

    # storage for controls discretized on intervals of tlist
    pulses0 :: Vector{Vector{Float64}}

    # second pulse storage: pulses0 and pulses1 alternate in storing the guess
    # pulses and optimized pulses in each iteration
    pulses1 ::  Vector{Vector{Float64}}

    # values of ∫gₐ(t)dt for each pulse
    g_a_int :: Vector{Float64}

    # update shapes S(t) for each pulse, discretized on intervals
    update_shapes ::  Vector{Vector{Float64}}

    lambda_vals :: Vector{Float64}

    # map of controls to options
    pulse_options :: POT  # TODO: this is not a good name

    # Result object

    result :: KrotovResult{STST}

    #################################
    # scratch objects, per objective:

    # backward-propagated states
    # note: storage for fw-propagated states is in result.states
    bw_states :: Vector{STST}

    # dynamical generator at a particular point in time
    G :: Vector{GT}

    vals_dict :: Vector{VDT}

    fw_storage :: Vector{STORT}  # forward storage array (per objective)

    fw_storage2 :: Vector{STORT}  # forward storage array (per objective)

    bw_storage :: Vector{STORT}  # backward storage array (per objective)

    prop_wrk :: Vector{PRWT}

    use_threads :: Bool

    function KrotovWrk(problem::QuantumControlBase.ControlProblem)
        prop_method = get(problem.kwargs, :prop_method, Val(:auto))
        use_threads = get(problem.kwargs, :use_threads, false)
        objectives = [obj for obj in problem.objectives]
        adjoint_objectives = [adjoint(obj) for obj in problem.objectives]
        controls = getcontrols(objectives)
        tlist = problem.tlist
        kwargs = problem.kwargs
        # TODO: check that kwargs has all the required parameters
        pulses0 = [
            discretize_on_midpoints(control, tlist) for control in controls
        ]
        pulses1 = [copy(pulse) for pulse in pulses0]
        g_a_int = zeros(length(pulses0))
        pulse_options = problem.pulse_options
        update_shapes = [
            discretize_on_midpoints(
                    pulse_options[control][:update_shape],
                    tlist
            )
            for control in controls
        ]
        lambda_vals = [
            pulse_options[control][:lambda_a] for control in controls
        ]
        result = KrotovResult(problem)
        bw_states = [similar(obj.initial_state) for obj in objectives]
        zero_vals = IdDict(control => zero(pulses0[i][1]) for (i, control) in enumerate(controls))
        G = [setcontrolvals(obj.generator, zero_vals) for obj in objectives]
        vals_dict = [copy(zero_vals) for _ in objectives]
        # TODO: forward_storage only if g_b != 0
        fw_storage = [init_storage(obj.initial_state, tlist) for obj in objectives]
        # TODO: second forward storage only if second order
        fw_storage2 = [init_storage(obj.initial_state, tlist) for obj in objectives]
        bw_storage = [init_storage(obj.initial_state, tlist) for obj in objectives]
        prop_wrk = [
            QuantumControlBase.initobjpropwrk(obj, tlist, prop_method; kwargs...)
            for obj in objectives
        ]
        # TODO: separate propwrk for backward propagation
        new{eltype(objectives),eltype(adjoint_objectives), typeof(kwargs),
            typeof(controls), typeof(pulse_options),
            typeof(objectives[1].initial_state), eltype(vals_dict),
            eltype(fw_storage), eltype(prop_wrk), eltype(G)
        }(
            objectives, adjoint_objectives, kwargs, controls, pulses0, pulses1,
            g_a_int, update_shapes, lambda_vals, pulse_options, result,
            bw_states, G, vals_dict, fw_storage, fw_storage2, bw_storage,
            prop_wrk, use_threads
        )
    end

end


"""Use Krotov's method to optimize the given optimization problem.

```julia
result = optimize_pulses(problem)
```

optimizes the control `problem`, see
[`QuantumControlBase.ControlProblem`](@ref).

Parameters are taken from the keyword arguments used in the instantiation of
`problem`.

# Required problem keyword arguments

The optimization functional is given implicitly via the mandatory `problem`
keyword argument `chi`.

# Optional problem keyword arguments

The following `problem` keyword arguments are supported (with default values):

* `sigma=nothing`: Function that calculate the second-order contribution. If
   not given, the first-order Krotov method is used.
* `iter_start=0`: the initial iteration number
* `iter_stop=5000`: the maximum iteration number
* `check_convergence`: a function to check whether convergence has been
  reached. Receives a [`KrotovResult`](@ref) object `result`, and should set
  `result.converged` to `true` and `result.message` to an appropriate string in
  case of convergence. Multiple convergence checks can be performed by chaining
  functions with `∘`.
* `prop_method=:auto`: The propagation method to use

"""
function optimize_pulses(problem)
    sigma = get(problem.kwargs, :sigma, nothing)
    iter_start = get(problem.kwargs, :iter_start, 0)
    info_hook = get(problem.kwargs, :info_hook, print_table)
    check_convergence! = get(problem.kwargs, :check_convergence, res -> res)
    # note: the default `check_convergence!` is a no-op. We still always check
    # for "Reached maximum number of iterations" in `update_result!`
    skip_initial_forward_propagation = get(
        problem.kwargs, :skip_initial_forward_propagation, false
    )

    wrk = KrotovWrk(problem)
    # TODO: if continuing previous optimization (`continue_from` argument),
    # ammend wrk from existing Result

    i = iter_start

    ϵ⁽ⁱ⁾ = wrk.pulses0
    ϵ⁽ⁱ⁺¹⁾ = wrk.pulses1

    if skip_initial_forward_propagation
        @info "Skipping initial forward propagation"
    else
        @threadsif wrk.use_threads for (k, obj) in collect(enumerate(wrk.objectives))
            krotov_initial_fw_prop!(
                ϵ⁽ⁱ⁾, wrk.result.states[k], obj.initial_state, k, wrk
            )
        end
    end

    # TODO: if sigma, fw_storage0 = fw_storage
    update_result!(wrk, 0)
    info_hook(wrk, 0)
    # TODO: update_hook

    while !wrk.result.converged
        i = i + 1
        krotov_iteration(wrk, ϵ⁽ⁱ⁾, ϵ⁽ⁱ⁺¹⁾)
        update_result!(wrk, i)
        info_hook(wrk, i)
        check_convergence!(wrk.result)
        # TODO: update_hook, e.g. to re-initialize propagation, apply spectral
        # filter
        ϵ⁽ⁱ⁾, ϵ⁽ⁱ⁺¹⁾ = ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾
    end

    finalize_result!(ϵ⁽ⁱ⁾, wrk)

    return wrk.result

end


function krotov_initial_fw_prop!(ϵ⁽⁰⁾, ϕₖ, ϕₖⁱⁿ, k, wrk)
    Φ₀ = wrk.fw_storage[k]
    copyto!(ϕₖ,  ϕₖⁱⁿ)
    (Φ₀ ≠ nothing) && write_to_storage!(Φ₀, 1,  ϕₖⁱⁿ)
    N_T = length(wrk.result.tlist) - 1
    for n = 1:N_T
        G, dt = _fw_gen(ϵ⁽⁰⁾, k, n, wrk)
        propstep!(ϕₖ, G, dt, wrk.prop_wrk[k])
        (Φ₀ ≠ nothing) && write_to_storage!(Φ₀, n+1, ϕₖ)
    end
    # TODO: allow a custom propstep! routine
end


function krotov_iteration(wrk, ϵ⁽ⁱ⁾, ϵ⁽ⁱ⁺¹⁾)

    ϕ = wrk.result.states  # assumed to contain the results of forward propagation
    χ = wrk.bw_states
    chi! = wrk.kwargs[:chi]
    N_T = length(wrk.result.tlist) - 1
    N = length(wrk.objectives)
    L = length(wrk.controls)
    X = wrk.bw_storage
    Φ = wrk.fw_storage  # TODO: pass in Φ₁, Φ₀ as parameters
    mu = get(wrk.kwargs, :mu, derivative_wrt_pulse)
    ∫gₐdt = wrk.g_a_int

    # TODO: re-initialize wrk.prop_wrk?

    # backward propagation
    chi!(χ, ϕ, wrk)
    @threadsif wrk.use_threads for k = 1:N
        write_to_storage!(X[k], N_T+1, χ[k])
        for n = N_T:-1:1
            local (G, dt) = _bw_gen(ϵ⁽ⁱ⁾, k, n, wrk)
            propstep!(χ[k], G, dt, wrk.prop_wrk[k])
            write_to_storage!(X[k], n, χ[k])
        end
    end

    # pulse update and forward propagation

    @threadsif wrk.use_threads for k = 1:N
        copyto!(ϕ[k], wrk.objectives[k].initial_state)
    end

    ∫gₐdt .= 0.0
    for n = 1:N_T  # `n` is the index for the time interval
        dt = wrk.result.tlist[n+1] - wrk.result.tlist[n]
        for k = 1:N
            get_from_storage!(χ[k], X[k], n)
        end
        Δuₙ = zeros(L)
        for l = 1:L  # `l` is the index for the different controls
            Sₗ = wrk.update_shapes[l]
            λₐ = wrk.lambda_vals[l]
            for k = 1:N
                μₗₖₙ = mu(wrk, l, k, n)
                Δuₙ[l] += (Sₗ[n]/λₐ) * imag(dot(χ[k], μₗₖₙ, ϕ[k]))
            end
            (∫gₐdt)[l] += abs(Δuₙ[l])^2 * dt
        end
        # TODO: second order update
        Δϵₙ = Δuₙ  # no parameterization (TODO)
        for l = 1:L
            ϵ⁽ⁱ⁺¹⁾[l][n] = ϵ⁽ⁱ⁾[l][n] + Δϵₙ[l]
        end
        @threadsif wrk.use_threads for k = 1:N
            local (G, dt) = _fw_gen(ϵ⁽ⁱ⁺¹⁾, k, n, wrk)
            propstep!(ϕ[k], G, dt, wrk.prop_wrk[k])
            write_to_storage!(Φ[k], n, ϕ[k])
        end
        # TODO: update sigma
    end
end


# The dynamical generator for the forward propagation
function _fw_gen(ϵ, k, n, wrk)
    vals_dict = wrk.vals_dict[k]
    t = wrk.result.tlist
    for (l, control) in enumerate(wrk.controls)
        vals_dict[control] = ϵ[l][n]
    end
    dt = t[n+1] - t[n]
    setcontrolvals!(wrk.G[k], wrk.objectives[k].generator, vals_dict)
    return wrk.G[k], dt
end


# The dynamical generator for the backward propagation
function _bw_gen(ϵ, k, n, wrk)
    vals_dict = wrk.vals_dict[k]
    t = wrk.result.tlist
    for (l, control) in enumerate(wrk.controls)
        vals_dict[control] = ϵ[l][n]
    end
    dt = t[n+1] - t[n]
    setcontrolvals!(wrk.G[k], wrk.adjoint_objectives[k].generator, vals_dict)
    return wrk.G[k], -dt
end


function derivative_wrt_pulse(wrk::KrotovWrk, l::Int64, k::Int64, n::Int64)
    # l: control index; k: objectives index, n: time grid interval index
    G = wrk.objectives[k].generator
    control = wrk.controls[l]
    return derivate_wrt_pulse(wrk, G, control, n)
end

function derivate_wrt_pulse(wrk::KrotovWrk, G::Tuple, control, n::Int64)
    # assume G to be a nested tuple, and that `control` only occurs once
    for part in G
        if isa(part, Tuple)
            if part[2] === control
                return part[1]
            end
        end
    end
end


function update_result!(wrk::KrotovWrk, i::Int64)
    res = wrk.result
    J_T_func = wrk.kwargs[:J_T]
    res.J_T_prev = res.J_T
    res.J_T = J_T_func(res.states, wrk)
    res.iter = i
    if i >= res.iter_stop
        res.converged = true
        res.message = "Reached maximum number of iterations"
        # Note: other convergence checks are done in user-supplied
        # check_convergence routine
    end
    prev_time = res.end_local_time
    res.end_local_time = now()
    res.secs = Dates.toms(res.end_local_time - prev_time) / 1000.0
    # TODO: calculate τ values
end


function finalize_result!(ϵ_opt, wrk::KrotovWrk)
    res = wrk.result
    res.end_local_time = now()
    for l in 1:length(ϵ_opt)
        res.optimized_controls[l] = discretize(ϵ_opt[l], res.tlist)
    end
end


"""Default info_hook"""
function print_table(wrk, iteration)
    J_T = wrk.result.J_T
    g_a_int = sum(wrk.g_a_int)
    J = J_T + g_a_int
    ΔJ_T = J_T - wrk.result.J_T_prev
    ΔJ = ΔJ_T + g_a_int
    secs = wrk.result.secs

    iter_stop = "$(get(wrk.kwargs, :iter_stop, 5000))"
    widths = [max(length("$iter_stop"), 6), 11, 11, 11, 11, 11, 8]

    if iteration == 0
        header = ["iter.", "J_T", "∫gₐ(t)dt", "J", "ΔJ_T", "ΔJ", "secs"]
        for (header, w) in zip(header, widths)
            print(lpad(header, w))
        end
        print("\n")
    end

    strs = (
        "$iteration",
        @sprintf("%.2e", J_T),
        @sprintf("%.2e", g_a_int),
        @sprintf("%.2e", J),
        (iteration > 0) ? @sprintf("%.2e", ΔJ_T) : "n/a",
        (iteration > 0) ? @sprintf("%.2e", ΔJ) : "n/a",
        @sprintf("%.1f", secs),
    )
    for (str, w) in zip(strs, widths)
        print(lpad(str, w))
    end
    print("\n")
end
