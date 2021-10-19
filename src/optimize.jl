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
        iter = iter_start
        secs = 0
        tau_vals = Vector{ComplexF64}()
        guess_controls = [
            discretize(control, tlist) for control in controls
        ]
        J_T = 0.0
        J_T_prev = 0.0
        optimized_controls = [copy(guess) for guess in guess_controls]
        states = [similar(obj.initial_state) for obj in problem.objectives]
        start_local_time = now()
        end_local_time = now()
        records = Vector{Tuple}()
        converged = false
        message = "in progress"
        new{eltype(states)}(
            tlist, iter_start, iter_stop, iter, secs, tau_vals, J_T, J_T_prev,
            guess_controls, optimized_controls, states, start_local_time,
            end_local_time, records, converged, message)
    end
end

Base.show(io::IO, r::KrotovResult) = print(io, "KrotovResult<$(r.message)>")
Base.show(io::IO, ::MIME"text/plain", r::KrotovResult) = print(io, """
Krotov Optimization Result
--------------------------
- Started at $(r.start_local_time)
- Number of objectives: $(length(r.states))
- Number of iterations: $(max(r.iter - r.iter_start, 0))
- Value of functional: $(r.J_T)
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

    is_parametrized :: Vector{Bool}

    parametrization :: Vector{PulseParametrization} # TODO

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

    control_derivs :: Vector{Vector{Union{Function, Nothing}}}

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
        control_derivs = [
            getcontrolderivs(obj.generator, controls) for obj in objectives
        ]
        tlist = problem.tlist
        kwargs = Dict(problem.kwargs)
        # TODO: check that kwargs has all the required parameters
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
        is_parametrized = [
            haskey(pulse_options[control], :parametrization)
            for control in controls
        ]
        parametrization = [
            get(pulse_options[control], :parametrization, NoParametrization())
            for control in controls
        ]
        if haskey(kwargs, :continue_from)
            @info "Continuing previous optimization"
            result = kwargs[:continue_from]
            if !(result isa KrotovResult)
                # account for continuing from a different optimization method
                result = convert(KrotovResult, result)
            end
            result.iter_stop = get(problem.kwargs, :iter_stop, 5000)
            result.converged = false
            result.start_local_time = now()
            result.message = "in progress"
            pulses0 = [
                discretize_on_midpoints(control, tlist)
                for control in result.optimized_controls
            ]
        else
            result = KrotovResult(problem)
            pulses0 = [
                discretize_on_midpoints(control, tlist) for control in controls
            ]
        end
        pulses1 = [copy(pulse) for pulse in pulses0]
        g_a_int = zeros(length(pulses0))
        bw_states = [similar(obj.initial_state) for obj in objectives]
        zero_vals = IdDict(control => zero(pulses0[i][1]) for (i, control) in enumerate(controls))
        G = [evalcontrols(obj.generator, zero_vals) for obj in objectives]
        vals_dict = [copy(zero_vals) for _ in objectives]
        # TODO: forward_storage only if g_b != 0
        fw_storage = [init_storage(obj.initial_state, tlist) for obj in objectives]
        # TODO: second forward storage only if second order
        fw_storage2 = [init_storage(obj.initial_state, tlist) for obj in objectives]
        bw_storage = [init_storage(obj.initial_state, tlist) for obj in objectives]
        prop_wrk = [
            QuantumControlBase.initobjpropwrk(obj, tlist, prop_method;
                                              initial_state=obj.initial_state,
                                              kwargs...)
            for obj in objectives
        ]
        # TODO: separate propwrk for backward propagation
        new{eltype(objectives),eltype(adjoint_objectives), typeof(kwargs),
            typeof(controls), typeof(pulse_options),
            typeof(objectives[1].initial_state), eltype(vals_dict),
            eltype(fw_storage), eltype(prop_wrk), eltype(G)
        }(
            objectives, adjoint_objectives, kwargs, controls, pulses0, pulses1,
            g_a_int, update_shapes, lambda_vals, is_parametrized,
            parametrization, pulse_options, result, bw_states, G,
            control_derivs, vals_dict, fw_storage, fw_storage2,
            bw_storage, prop_wrk, use_threads
        )
    end

end


"""Use Krotov's method to optimize the given optimization problem.

```julia
result = optimize_pulses(problem; kwargs...)
```

optimizes the given control problem, see
[`QuantumControlBase.ControlProblem`](@ref).

Keyword arguments that control the optimization are taken from the keyword
arguments used in the instantiation of `problem`. Any `kwargs` passed directly
to `optimize_pulses` will update (overwrite) the parameters in `problem`.

# Required problem keyword arguments

* `J_T`: A function `J_T(ϕ, objectives)` that evaluates the final time
  functional from a list `ϕ` of forward-propagated states and
  `problem.objectives`.
* `chi`: A function `chi!(χ, ϕ, objectives)` what receives a list `ϕ`
  of the forward propagates state and must set ``χₖ=∂J_T/∂⟨ϕₖ|``.

# Optional problem keyword arguments

The following keyword arguments are supported (with default values):

* `sigma=nothing`: Function that calculate the second-order contribution. If
   not given, the first-order Krotov method is used.
* `iter_start=0`: the initial iteration number
* `iter_stop=5000`: the maximum iteration number
* `prop_method=:auto`: The propagation method to use
* `update_hook`: A function that receives the Krotov workspace, the iteration
   number, the list of updated pulses and the list of guess pulses as
   positional arguments. The function may mutate any of its arguments. This may
   be used e.g. to apply a spectral filter to the updated pulses, or to update
   propagation workspaces inside the Krotov workspace.
* `info_hook`: A function that receives the same argumens as `update_hook`, in
   order to write information about the current iteration to the screen or to a
   file. The default `info_hook` prints a table with convergence information to
   the screen. Runs after `update_hook`. The `info_hook` function may return a
   tuple, which is stored in the list of `records` inside the
   [`KrotovResult`](@ref) object.
* `check_convergence`: a function to check whether convergence has been
  reached. Receives a [`KrotovResult`](@ref) object `result`, and should set
  `result.converged` to `true` and `result.message` to an appropriate string in
  case of convergence. Multiple convergence checks can be performed by chaining
  functions with `∘`. The convergence check is performed after any calls to
  `update_hook` and `info_hook`.

"""
function optimize_pulses(problem; kwargs...)
    merge!(problem.kwargs, kwargs)
    sigma = get(problem.kwargs, :sigma, nothing)
    iter_start = get(problem.kwargs, :iter_start, 0)
    update_hook! = get(problem.kwargs, :update_hook, (args...) -> nothing)
    info_hook = get(problem.kwargs, :info_hook, print_table)
    check_convergence! = get(problem.kwargs, :check_convergence, res -> res)
    # note: the default `check_convergence!` is a no-op. We still always check
    # for "Reached maximum number of iterations" in `update_result!`
    skip_initial_forward_propagation = get(
        problem.kwargs, :skip_initial_forward_propagation, false
    )

    wrk = KrotovWrk(problem)

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
    update_hook!(wrk, 0, ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾)
    info_tuple = info_hook(wrk, 0, ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾)
    (info_tuple !== nothing) && push!(wrk.records, info_tuple)

    i = wrk.result.iter  # = 0, unless continuing from previous optimization
    while !wrk.result.converged
        i = i + 1
        krotov_iteration(wrk, ϵ⁽ⁱ⁾, ϵ⁽ⁱ⁺¹⁾)
        update_result!(wrk, i)
        update_hook!(wrk, i, ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾)
        info_tuple = info_hook(wrk, i, ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾)
        (info_tuple !== nothing) && push!(wrk.records, info_tuple)
        check_convergence!(wrk.result)
        ϵ⁽ⁱ⁾, ϵ⁽ⁱ⁺¹⁾ = ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾
    end

    finalize_result!(ϵ⁽ⁱ⁾, wrk)

    return wrk.result

end


function krotov_initial_fw_prop!(ϵ⁽⁰⁾, ϕₖ, ϕₖⁱⁿ, k, wrk)
    Φ₀ = wrk.fw_storage[k]
    copyto!(ϕₖ,  ϕₖⁱⁿ)
    (Φ₀ !== nothing) && write_to_storage!(Φ₀, 1,  ϕₖⁱⁿ)
    N_T = length(wrk.result.tlist) - 1
    for n = 1:N_T
        G, dt = _fw_gen(ϵ⁽⁰⁾, k, n, wrk)
        propstep!(ϕₖ, G, dt, wrk.prop_wrk[k])
        (Φ₀ !== nothing) && write_to_storage!(Φ₀, n+1, ϕₖ)
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
    ∫gₐdt = wrk.g_a_int
    Im = imag
    # TODO: re-initialize wrk.prop_wrk?

    # backward propagation
    chi!(χ, ϕ, wrk.objectives)
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
        ϵₙ⁽ⁱ⁺¹⁾ = [ϵ⁽ⁱ⁾[l][n]  for l ∈ 1:L]  # ϵₙ⁽ⁱ⁺¹⁾ ≈ ϵₙ⁽ⁱ⁾ for non-linear controls
        # TODO: we could add a self-consistent loop here for ϵₙ⁽ⁱ⁺¹⁾
        Δuₙ = zeros(L)
        for l = 1:L  # `l` is the index for the different controls
            Sₗ = wrk.update_shapes[l]
            λₐ = wrk.lambda_vals[l]
            ∂ϵₗ╱∂u :: Function = wrk.parametrization[l].de_du_derivative
            uₗ = wrk.parametrization[l].u_of_epsilon
            for k = 1:N  # k is the index over the objectives
                ∂Hₖ╱∂ϵₗ :: Union{Function, Nothing} = wrk.control_derivs[k][l]
                if !isnothing(∂Hₖ╱∂ϵₗ)
                    μₗₖₙ = (∂Hₖ╱∂ϵₗ)(ϵₙ⁽ⁱ⁺¹⁾[l])
                    αₗ = (Sₗ[n]/λₐ)  # Krotov step size
                    if wrk.is_parametrized[l]
                        ∂ϵₗ╱∂uₗ = (∂ϵₗ╱∂u)(uₗ(ϵₙ⁽ⁱ⁺¹⁾[l]))
                        Δuₙ[l] += αₗ * ∂ϵₗ╱∂uₗ * Im(dot(χ[k], μₗₖₙ, ϕ[k]))
                    else
                        Δuₙ[l] += αₗ * Im(dot(χ[k], μₗₖₙ, ϕ[k]))
                    end
                end
            end
        end
        # TODO: second order update
        for l = 1:L
            if wrk.is_parametrized[l]
                uₗ = wrk.parametrization[l].u_of_epsilon
                ϵₗ = wrk.parametrization[l].epsilon_of_u
                ϵ⁽ⁱ⁺¹⁾[l][n] = ϵₗ(uₗ(ϵ⁽ⁱ⁾[l][n]) + Δuₙ[l])
            else
                Δϵₗₙ = Δuₙ[l]
                ϵ⁽ⁱ⁺¹⁾[l][n] = ϵ⁽ⁱ⁾[l][n] + Δϵₗₙ
            end
        end
        # TODO: end of self-consistent loop
        @. ∫gₐdt += abs(Δuₙ)^2 * dt
        @threadsif wrk.use_threads for k = 1:N
            local (G, dt) = _fw_gen(ϵ⁽ⁱ⁺¹⁾, k, n, wrk)
            propstep!(ϕ[k], G, dt, wrk.prop_wrk[k])
            write_to_storage!(Φ[k], n, ϕ[k])
        end
        # TODO: update sigma
    end  # time loop
end


# The dynamical generator for the forward propagation
function _fw_gen(ϵ, k, n, wrk)
    vals_dict = wrk.vals_dict[k]
    t = wrk.result.tlist
    for (l, control) in enumerate(wrk.controls)
        vals_dict[control] = ϵ[l][n]
    end
    dt = t[n+1] - t[n]
    evalcontrols!(wrk.G[k], wrk.objectives[k].generator, vals_dict)
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
    evalcontrols!(wrk.G[k], wrk.adjoint_objectives[k].generator, vals_dict)
    return wrk.G[k], -dt
end


function update_result!(wrk::KrotovWrk, i::Int64)
    res = wrk.result
    J_T_func = wrk.kwargs[:J_T]
    res.J_T_prev = res.J_T
    res.J_T = J_T_func(res.states, wrk.objectives)
    (i > 0) && (res.iter = i)
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


"""Print optimization progress as a table.

This functions serves as the default `info_hook` for an optimization with
Krotov's method.
"""
function print_table(wrk, iteration, args...)
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
