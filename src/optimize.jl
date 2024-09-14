using QuantumControl.QuantumPropagators.Generators: Operator
using QuantumControl.QuantumPropagators.Controls: discretize, evaluate
using QuantumControl.QuantumPropagators.Interfaces: supports_inplace
using QuantumControl.QuantumPropagators: prop_step!, reinit_prop!, propagate
using QuantumControl.QuantumPropagators.Storage:
    write_to_storage!, get_from_storage!, get_from_storage
using QuantumControl.Functionals: make_chi, taus!
using QuantumControl: set_atexit_save_optimization
using QuantumControl: @threadsif, Trajectory
using QuantumControl.QuantumPropagators: _StoreState
using LinearAlgebra
using Printf

import QuantumControl: optimize, make_print_iters

@doc raw"""
```julia
using Krotov
result = optimize(problem; method=Krotov, kwargs...)
```

optimizes the given control [`problem`](@ref QuantumControl.ControlProblem)
using Krotov's method, returning a [`KrotovResult`](@ref).

Keyword arguments that control the optimization are taken from the keyword
arguments used in the instantiation of `problem`; any of these can be overridden
with explicit keyword arguments to `optimize`.

# Required problem keyword arguments

* `J_T`: A function `J_T(Ψ, trajectories)` that evaluates the final time
  functional from a list `Ψ` of forward-propagated states and
  `problem.trajectories`. The function `J_T` may also take a keyword argument
  `tau`. If it does, a vector containing the complex overlaps of the target
  states (`target_state` property of each trajectory in `problem.trajectories`)
  with the propagated states will be passed to `J_T`.

# Recommended problem keyword arguments

* `lambda_a=1.0`: The inverse Krotov step width λₐ for every pulse.
* `update_shape=(t->1.0)`: A function `S(t)` for the "update shape" that scales
  the update for every pulse.

If different controls require different `lambda_a` or `update_shape`, a dict
`pulse_options` must be given instead of a global `lambda_a` and
`update_shape`; see below.

# Optional problem keyword arguments

The following keyword arguments are supported (with default values):

* `pulse_options`: A dictionary that maps every control (as obtained by
  [`get_controls`](@ref
  QuantumControl.QuantumPropagators.Controls.get_controls) from the
  `problem.trajectories`) to the following dict:

  - `:lambda_a`:  The value for inverse Krotov step width λₐ.
  - `:update_shape`: A function `S(t)` for the "update shape" that scales
    the Krotov pulse update.

  This overrides the global `lambda_a` and `update_shape` arguments.

* `chi`: A function `chi(Ψ, trajectories)` that receives a list `Ψ`
  of the forward propagated states and returns a vector of states
  ``|χₖ⟩ = -∂J_T/∂⟨Ψₖ|``. If not given, it will be automatically determined
  from `J_T` via [`make_chi`](@ref) with the default parameters. Similarly to
  `J_T`, if `chi` accepts a keyword argument `tau`, it will be passed a vector
  of complex overlaps.
* `sigma=nothing`: A function that calculates the second-order contribution. If
  not given, the first-order Krotov method is used.
* `iter_start=0`: The initial iteration number.
* `iter_stop=5000`: The maximum iteration number.
* `prop_method`: The propagation method to use for each trajectory; see below.
* `print_iters=true`: Whether to print information after each iteration.
* `store_iter_info=Set()`: Which fields from `print_iters` to store in
  `result.records`. A subset of
  `Set(["iter.", "J_T", "∫gₐ(t)dt", "J", "ΔJ_T", "ΔJ", "secs"])`.
* `callback`: A function (or tuple of functions) that receives the
  [Krotov workspace](@ref KrotovWrk), the iteration number, the list of updated
  pulses, and the list of guess pulses as positional arguments. The function
  may return a tuple of values which are stored in the
  [`KrotovResult`](@ref) object `result.records`. The function can also mutate
  any of its arguments, in particular the updated pulses. This may be used,
  e.g., to apply a spectral filter to the updated pulses or to perform
  similar manipulations. Note that `print_iters=true` (default) adds an
  automatic callback to print information after each iteration. With
  `store_iter_info`, that callback automatically stores a subset of the
  printed information.
* `check_convergence`: A function to check whether convergence has been
  reached. Receives a [`KrotovResult`](@ref) object `result`, and should set
  `result.converged` to `true` and `result.message` to an appropriate string in
  case of convergence. Multiple convergence checks can be performed by chaining
  functions with `∘`. The convergence check is performed after any `callback`.
* `verbose=false`: If `true`, print information during initialization.
* `rethrow_exceptions`: By default, any exception ends the optimization but
  still returns a [`KrotovResult`](@ref) that captures the message associated
  with the exception. This is to avoid losing results from a long-running
  optimization when an exception occurs in a later iteration. If
  `rethrow_exceptions=true`, instead of capturing the exception, it will be
  thrown normally.

# Trajectory propagation

Krotov's method involves the forward and backward propagation for every
[`Trajectory`](@ref) in the `problem`. The keyword arguments for each
propagation (see [`propagate`](@ref)) are determined from any properties of
each [`Trajectory`](@ref) that have a `prop_` prefix, cf.
[`init_prop_trajectory`](@ref).

In situations where different parameters are required for the forward and
backward propagation, instead of the `prop_` prefix, the `fw_prop_` and
`bw_prop_` prefixes can be used, respectively. These override any setting with
the `prop_` prefix. This applies both to the properties of each
[`Trajectory`](@ref) and the problem keyword arguments.

Note that the propagation method for each propagation must be specified. In
most cases, it is sufficient (and recommended) to pass a global `prop_method`
problem keyword argument.
"""
optimize(problem, method::Val{:Krotov}) = optimize_krotov(problem)
optimize(problem, method::Val{:krotov}) = optimize_krotov(problem)

"""
See [`optimize(problem; method=Krotov, kwargs...)`](@ref optimize(::Any, ::Val{:krotov})).
"""
function optimize_krotov(problem)
    callback = get(problem.kwargs, :callback, (args...) -> nothing)
    if haskey(problem.kwargs, :update_hook) || haskey(problem.kwargs, :info_hook)
        msg = "The `update_hook` and `info_hook` arguments have been superseded by the `callback` argument"
        throw(ArgumentError(msg))
    end
    check_convergence! = get(problem.kwargs, :check_convergence, res -> res)
    # note: the default `check_convergence!` is a no-op. We still always check
    # for "Reached maximum number of iterations" in `update_result!`
    verbose = get(problem.kwargs, :verbose, false)
    skip_initial_forward_propagation =
        get(problem.kwargs, :skip_initial_forward_propagation, false)

    wrk = KrotovWrk(problem; verbose)

    ϵ⁽ⁱ⁾ = wrk.pulses0
    ϵ⁽ⁱ⁺¹⁾ = wrk.pulses1

    if skip_initial_forward_propagation
        @info "Skipping initial forward propagation"
    else
        @threadsif wrk.use_threads for (k, traj) in collect(enumerate(wrk.trajectories))
            krotov_initial_fw_prop!(ϵ⁽ⁱ⁾, traj.initial_state, k, wrk)
        end
    end

    # TODO: if sigma, fw_storage0 = fw_storage
    update_result!(wrk, 0)
    info_tuple = callback(wrk, 0, ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾)
    if !(isnothing(info_tuple) || isempty(info_tuple))
        push!(wrk.result.records, info_tuple)
    end

    i = wrk.result.iter  # = 0, unless continuing from previous optimization
    atexit_filename = get(problem.kwargs, :atexit_filename, nothing)
    # atexit_filename is undocumented on purpose: this is considered a feature
    # of @optimize_or_load
    if !isnothing(atexit_filename)
        set_atexit_save_optimization(atexit_filename, wrk.result)
        if !isinteractive()
            @info "Set callback to store result in $(relpath(atexit_filename)) on unexpected exit."
            # In interactive mode, `atexit` is very unlikely, and
            # `InterruptException` is handles via try/catch instead.
        end
    end
    try
        while !wrk.result.converged
            i = i + 1
            krotov_iteration(wrk, ϵ⁽ⁱ⁾, ϵ⁽ⁱ⁺¹⁾)
            update_result!(wrk, i)
            info_tuple = callback(wrk, i, ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾)
            if !(isnothing(info_tuple) || isempty(info_tuple))
                push!(wrk.result.records, info_tuple)
            end
            check_convergence!(wrk.result)
            ϵ⁽ⁱ⁾, ϵ⁽ⁱ⁺¹⁾ = ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾
        end
    catch exc
        if get(problem.kwargs, :rethrow_exceptions, false)
            rethrow()
        end
        # Primarily, this is intended to catch Ctrl-C in interactive
        # optimizations (InterruptException)
        exc_msg = sprint(showerror, exc)
        wrk.result.message = "Exception: $exc_msg"
    end

    finalize_result!(ϵ⁽ⁱ⁾, wrk)
    if !isnothing(atexit_filename)
        popfirst!(Base.atexit_hooks)
    end

    return wrk.result

end


function transform_control_ranges(c, ϵ_min, ϵ_max, check)
    if check
        return (min(ϵ_min, 2 * ϵ_min), max(ϵ_max, 2 * ϵ_max))
    else
        return (min(ϵ_min, 5 * ϵ_min), max(ϵ_max, 5 * ϵ_max))
    end
end


function krotov_initial_fw_prop!(ϵ⁽⁰⁾, ϕₖ, k, wrk)
    for propagator in wrk.fw_propagators
        propagator.parameters = IdDict(zip(wrk.controls, ϵ⁽⁰⁾))
    end
    reinit_prop!(wrk.fw_propagators[k], ϕₖ; transform_control_ranges)

    Φ₀ = wrk.fw_storage[k]
    (Φ₀ !== nothing) && write_to_storage!(Φ₀, 1, ϕₖ)
    N_T = length(wrk.result.tlist) - 1
    for n = 1:N_T
        Ψₖ = prop_step!(wrk.fw_propagators[k])
        if haskey(wrk.fw_prop_kwargs[k], :callback)
            local cb = wrk.fw_prop_kwargs[k][:callback]
            local observables = get(wrk.fw_prop_kwargs[k], :observables, _StoreState())
            cb(wrk.fw_propagators[k], observables)
        end
        (Φ₀ !== nothing) && write_to_storage!(Φ₀, n + 1, Ψₖ)
    end
end


function _eval_mu(μ, wrk, ϵₙ, tlist, n)
    # Implementation for non-linear control terms
    vals_dict = IdDict(control => val for (control, val) ∈ zip(wrk.controls, ϵₙ))
    evaluate(μ, tlist, n; vals_dict)
end

# Linear control terms (static derivatives)
_eval_mu(μ::Operator, _...) = μ
_eval_mu(μ::AbstractMatrix, _...) = μ


function krotov_iteration(wrk, ϵ⁽ⁱ⁾, ϵ⁽ⁱ⁺¹⁾)

    tlist = wrk.result.tlist
    N_T = length(tlist) - 1
    N = length(wrk.trajectories)
    L = length(wrk.controls)
    X = wrk.bw_storage
    Φ = wrk.fw_storage  # TODO: pass in Φ₁, Φ₀ as parameters
    ∫gₐdt = wrk.g_a_int
    Im = imag
    chi = wrk.kwargs[:chi]  # guaranteed to exist in `KrotovWrk` constructor


    guess_parameters = IdDict(zip(wrk.controls, ϵ⁽ⁱ⁾))
    updated_parameters = IdDict(zip(wrk.controls, ϵ⁽ⁱ⁺¹⁾))

    # backward propagation

    Ψ = [propagator.state for propagator in wrk.fw_propagators]
    if wrk.chi_takes_tau
        χ = chi(Ψ, wrk.trajectories; tau=wrk.result.tau_vals)
    else
        χ = chi(Ψ, wrk.trajectories)
    end
    @threadsif wrk.use_threads for k = 1:N
        # TODO: normalize χ; warn if norm is close to zero
        wrk.bw_propagators[k].parameters = guess_parameters
        reinit_prop!(wrk.bw_propagators[k], χ[k]; transform_control_ranges)
        write_to_storage!(X[k], N_T + 1, χ[k])
        for n = N_T:-1:1
            local χₖ = prop_step!(wrk.bw_propagators[k])
            if haskey(wrk.bw_prop_kwargs[k], :callback)
                local cb = wrk.bw_prop_kwargs[k][:callback]
                local observables = get(wrk.bw_prop_kwargs[k], :observables, _StoreState())
                cb(wrk.bw_propagators[k], observables)
            end
            write_to_storage!(X[k], n, χₖ)
        end
    end

    # pulse update and forward propagation

    @threadsif wrk.use_threads for k = 1:N
        wrk.fw_propagators[k].parameters = updated_parameters
        local Ψₖ = wrk.trajectories[k].initial_state
        reinit_prop!(wrk.fw_propagators[k], Ψₖ; transform_control_ranges)
    end

    ∫gₐdt .= 0.0
    for n = 1:N_T  # `n` is the index for the time interval
        dt = tlist[n+1] - tlist[n]
        for k = 1:N
            if supports_inplace(χ[k])
                get_from_storage!(χ[k], X[k], n)
            else
                χ[k] = get_from_storage(X[k], n)
            end
        end
        ϵₙ⁽ⁱ⁺¹⁾ = [ϵ⁽ⁱ⁾[l][n] for l ∈ 1:L]  # ϵₙ⁽ⁱ⁺¹⁾ ≈ ϵₙ⁽ⁱ⁾ for non-linear controls
        # TODO: we could add a self-consistent loop here for ϵₙ⁽ⁱ⁺¹⁾
        Δuₙ = zeros(L)  # Δu is Δϵ without (Sₗₙ/λₐ) factor
        for l = 1:L  # `l` is the index for the different controls
            for k = 1:N  # k is the index over the trajectories
                Ψₖ = wrk.fw_propagators[k].state
                μₖₗ = wrk.control_derivs[k][l]
                if !isnothing(μₖₗ)
                    μₗₖₙ = _eval_mu(μₖₗ, wrk, ϵₙ⁽ⁱ⁺¹⁾, tlist, n)
                    Δuₙ[l] += Im(dot(χ[k], μₗₖₙ, Ψₖ))
                end
            end
        end
        # TODO: second order update
        for l = 1:L
            Sₗ = wrk.update_shapes[l]
            λₐ = wrk.lambda_vals[l]
            αₗ = (Sₗ[n] / λₐ)  # Krotov step size
            Δϵₗₙ = αₗ * Δuₙ[l]
            ϵ⁽ⁱ⁺¹⁾[l][n] = ϵ⁽ⁱ⁾[l][n] + Δϵₗₙ
            (∫gₐdt)[l] += αₗ * abs(Δuₙ[l])^2 * dt
        end
        # TODO: end of self-consistent loop
        @threadsif wrk.use_threads for k = 1:N
            local Ψₖ = prop_step!(wrk.fw_propagators[k])
            if haskey(wrk.fw_prop_kwargs[k], :callback)
                local cb = wrk.fw_prop_kwargs[k][:callback]
                local observables = get(wrk.fw_prop_kwargs[k], :observables, _StoreState())
                cb(wrk.fw_propagators[k], observables)
            end
            write_to_storage!(Φ[k], n, Ψₖ)
        end
        # TODO: update sigma
    end  # time loop
end


function update_result!(wrk::KrotovWrk, i::Int64)
    res = wrk.result
    J_T = wrk.kwargs[:J_T]
    res.J_T_prev = res.J_T
    for k in eachindex(wrk.fw_propagators)
        res.states[k] = wrk.fw_propagators[k].state
    end
    taus!(res.tau_vals, res.states, wrk.trajectories; ignore_missing_target_state=true)
    if wrk.J_T_takes_tau
        res.J_T = J_T(res.states, wrk.trajectories; tau=res.tau_vals)
    else
        res.J_T = J_T(res.states, wrk.trajectories)
    end
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
end


function finalize_result!(ϵ_opt, wrk::KrotovWrk)
    res = wrk.result
    res.end_local_time = now()
    for l in eachindex(ϵ_opt)
        res.optimized_controls[l] = discretize(ϵ_opt[l], res.tlist)
    end
end


make_print_iters(::Val{:Krotov}; kwargs...) = make_krotov_print_iters(; kwargs...)
make_print_iters(::Val{:krotov}; kwargs...) = make_krotov_print_iters(; kwargs...)


function make_krotov_print_iters(; kwargs...)

    header = ["iter.", "J_T", "∫gₐ(t)dt", "J", "ΔJ_T", "ΔJ", "secs"]
    store_iter_info = Set(get(kwargs, :store_iter_info, Set()))
    info_vals = Vector{Any}(undef, length(header))
    fill!(info_vals, nothing)
    store_iter = false
    store_J_T = false
    store_g_a_int = false
    store_J = false
    store_ΔJ_T = false
    store_ΔJ = false
    store_secs = false
    for item in store_iter_info
        if item == "iter."
            store_iter = true
        elseif item == "J_T"
            store_J_T = true
        elseif item == "∫gₐ(t)dt"
            store_g_a_int = true
        elseif item == "J"
            store_J = true
        elseif item == "ΔJ_T"
            store_ΔJ_T = true
        elseif item == "ΔJ"
            store_ΔJ = true
        elseif item == "secs"
            store_secs = true
        else
            msg = "Item $(repr(item)) in `store_iter_info` is not one of $(repr(header)))"
            throw(ArgumentError(msg))
        end
    end


    function print_table(wrk, iteration, args...)

        J_T = wrk.result.J_T
        g_a_int = sum(wrk.g_a_int)
        J = J_T + g_a_int
        ΔJ_T = J_T - wrk.result.J_T_prev
        ΔJ = ΔJ_T + g_a_int
        secs = wrk.result.secs

        store_iter && (info_vals[1] = iteration)
        store_J_T && (info_vals[2] = J_T)
        store_g_a_int && (info_vals[3] = g_a_int)
        store_J && (info_vals[4] = J)
        store_ΔJ_T && (info_vals[5] = ΔJ_T)
        store_ΔJ && (info_vals[6] = ΔJ)
        store_secs && (info_vals[7] = secs)

        iter_stop = "$(get(wrk.kwargs, :iter_stop, 5000))"
        widths = [max(length("$iter_stop"), 6), 11, 11, 11, 11, 11, 8]

        if iteration == 0
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
        flush(stdout)

        return Tuple((value for value in info_vals if (value !== nothing)))

    end

    return print_table

end
