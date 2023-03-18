using QuantumControlBase.QuantumPropagators.Generators: Operator
using QuantumControlBase.QuantumPropagators.Controls: discretize, evaluate
using QuantumControlBase.QuantumPropagators: prop_step!, reinit_prop!
using QuantumControlBase.QuantumPropagators.Storage: write_to_storage!, get_from_storage!
using QuantumControlBase: make_chi, set_atexit_save_optimization
using QuantumControlBase: @threadsif
using LinearAlgebra
using Printf

import QuantumControlBase: optimize

@doc raw"""
```julia
result = optimize(problem; method=:krotov, kwargs...)
```

optimizes the given
control [`problem`](@ref QuantumControlBase.ControlProblem) using Krotov's
method, returning a [`KrotovResult`](@ref).

Keyword arguments that control the optimization are taken from the keyword
arguments used in the instantiation of `problem`.

# Required problem keyword arguments

* `J_T`: A function `J_T(ϕ, objectives)` that evaluates the final time
  functional from a list `ϕ` of forward-propagated states and
  `problem.objectives`.

# Recommended problem keyword arguments

* `lambda_a=1.0`: The inverse Krotov step width λ_a for every pulse.
* `update_shape=(t->1.0)`: A function `S(t)` for the "update shape" that scales
  the update for every pulse

If different controls require different `lambda_a` or `update_shape`, a dict
`pulse_options` must be given instead of a global `lambda_a` and
`update_shape`, see below.

# Optional problem keyword arguments

The following keyword arguments are supported (with default values):

* `pulse_options`: A dictionary that maps every control (as obtained by
  [`get_controls`](@ref
  QuantumControlBase.QuantumPropagators.Controls.get_controls) from the
  `problem.objectives`) to the following dict:

  - `:lambda_a`:  The value for inverse Krotov step width λₐ
  - `:update_shape`: A function `S(t)` for the "update shape" that scales
    the Krotov pulse update.

  This overrides the global `lambda_a` and `update_shape` arguments.

* `chi`: A function `chi!(χ, ϕ, objectives)` what receives a list `ϕ`
  of the forward propagated states and must set ``|χₖ⟩ = -∂J_T/∂⟨ϕₖ|``. If not
  given, it will be automatically determined from `J_T` via [`make_chi`](@ref
  QuantumControlBase.make_chi) with the default parameters.
* `sigma=nothing`: Function that calculate the second-order contribution. If
  not given, the first-order Krotov method is used.
* `iter_start=0`: the initial iteration number
* `iter_stop=5000`: the maximum iteration number
* `prop_method`/`fw_prop_method`/`bw_prop_method`: The propagation method to
  use for each objective, see below.
* `update_hook`: A function that receives the Krotov workspace, the iteration
  number, the list of updated pulses and the list of guess pulses as
  positional arguments. The function may mutate any of its arguments. This may
  be used e.g. to apply a spectral filter to the updated pulses, or to update
  propagation workspaces inside the Krotov workspace.
* `info_hook`: A function (or tuple of functions) that receives the same
  argumens as `update_hook`, in order to write information about the current
  iteration to the screen or to a file. The default `info_hook` prints a table
  with convergence information to the screen. Runs after `update_hook`. The
  `info_hook` function may return a tuple, which is stored in the list of
  `records` inside the [`KrotovResult`](@ref) object.
* `check_convergence`: a function to check whether convergence has been
  reached. Receives a [`KrotovResult`](@ref) object `result`, and should set
  `result.converged` to `true` and `result.message` to an appropriate string in
  case of convergence. Multiple convergence checks can be performed by chaining
  functions with `∘`. The convergence check is performed after any calls to
  `update_hook` and `info_hook`.
* `verbose=false`: If `true`, print information during initialization

The propagation method for the forward propagation of each objective is
determined by the first available item of the following:

* a `fw_prop_method` keyword argument
* a `prop_method` keyword argument
* a property `fw_prop_method` of the objective
* a property `prop_method` of the objective
* the value `:auto`

The propagation method for the backword propagation is determined similarly,
but with `bw_prop_method` instead of `fw_prop_method`.
"""
optimize(problem, method::Val{:krotov}) = optimize_krotov(problem)
optimize(problem, method::Val{:Krotov}) = optimize_krotov(problem)

"""
See [`optimize(problem; method=:krotov, kwargs...)`](@ref optimize(::Any, ::Val{:krotov})).
"""
function optimize_krotov(problem)
    sigma = get(problem.kwargs, :sigma, nothing)
    iter_start = get(problem.kwargs, :iter_start, 0)
    update_hook! = get(problem.kwargs, :update_hook, (args...) -> nothing)
    info_hook = get(problem.kwargs, :info_hook, print_table)
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
        @threadsif wrk.use_threads for (k, obj) in collect(enumerate(wrk.objectives))
            krotov_initial_fw_prop!(ϵ⁽ⁱ⁾, obj.initial_state, k, wrk)
        end
    end

    # TODO: if sigma, fw_storage0 = fw_storage
    update_result!(wrk, 0)
    update_hook!(wrk, 0, ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾)
    info_tuple = info_hook(wrk, 0, ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾)
    (info_tuple !== nothing) && push!(wrk.result.records, info_tuple)

    i = wrk.result.iter  # = 0, unless continuing from previous optimization
    atexit_filename = get(problem.kwargs, :atexit_filename, nothing)
    # atexit_filename is undocumented on purpose: this is considered a feature
    # of @optimize_or_load
    if !isnothing(atexit_filename)
        set_atexit_save_optimization(atexit_filename, wrk.result)
        if !isinteractive()
            @info "Set callback to store result in $atexit_filename on unexpected exit."
            # In interactive mode, `atexit` is very unlikely, and
            # `InterruptException` is handles via try/catch instead.
        end
    end
    try
        while !wrk.result.converged
            i = i + 1
            krotov_iteration(wrk, ϵ⁽ⁱ⁾, ϵ⁽ⁱ⁺¹⁾)
            update_result!(wrk, i)
            update_hook!(wrk, i, ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾)
            info_tuple = info_hook(wrk, i, ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾)
            (info_tuple !== nothing) && push!(wrk.result.records, info_tuple)
            check_convergence!(wrk.result)
            ϵ⁽ⁱ⁾, ϵ⁽ⁱ⁺¹⁾ = ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾
        end
    catch exc
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


function krotov_initial_fw_prop!(ϵ⁽⁰⁾, ϕₖⁱⁿ, k, wrk)
    for propagator in wrk.fw_propagators
        propagator.parameters = IdDict(zip(wrk.controls, ϵ⁽⁰⁾))
    end
    reinit_prop!(wrk.fw_propagators[k], ϕₖⁱⁿ; transform_control_ranges)

    Φ₀ = wrk.fw_storage[k]
    (Φ₀ !== nothing) && write_to_storage!(Φ₀, 1, ϕₖⁱⁿ)
    N_T = length(wrk.result.tlist) - 1
    for n = 1:N_T
        ϕₖ = prop_step!(wrk.fw_propagators[k])
        (Φ₀ !== nothing) && write_to_storage!(Φ₀, n + 1, ϕₖ)
    end
    # TODO: allow a custom prop_step! routine
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

    χ = [propagator.state for propagator in wrk.bw_propagators]
    J_T_func = wrk.kwargs[:J_T]
    chi! = get(wrk.kwargs, :chi, make_chi(J_T_func, wrk.objectives))
    tlist = wrk.result.tlist
    N_T = length(tlist) - 1
    N = length(wrk.objectives)
    L = length(wrk.controls)
    X = wrk.bw_storage
    Φ = wrk.fw_storage  # TODO: pass in Φ₁, Φ₀ as parameters
    ∫gₐdt = wrk.g_a_int
    Im = imag

    guess_parameters = IdDict(zip(wrk.controls, ϵ⁽ⁱ⁾))
    updated_parameters = IdDict(zip(wrk.controls, ϵ⁽ⁱ⁺¹⁾))

    # backward propagation
    ϕ = [propagator.state for propagator in wrk.fw_propagators]
    chi!(χ, ϕ, wrk.objectives)
    @threadsif wrk.use_threads for k = 1:N
        wrk.bw_propagators[k].parameters = guess_parameters
        reinit_prop!(wrk.bw_propagators[k], χ[k]; transform_control_ranges)
        write_to_storage!(X[k], N_T + 1, χ[k])
        for n = N_T:-1:1
            local χₖ = prop_step!(wrk.bw_propagators[k])
            write_to_storage!(X[k], n, χₖ)
        end
    end

    # pulse update and forward propagation

    @threadsif wrk.use_threads for k = 1:N
        wrk.fw_propagators[k].parameters = updated_parameters
        local ϕₖ = wrk.objectives[k].initial_state
        reinit_prop!(wrk.fw_propagators[k], ϕₖ; transform_control_ranges)
    end

    ∫gₐdt .= 0.0
    for n = 1:N_T  # `n` is the index for the time interval
        dt = tlist[n+1] - tlist[n]
        for k = 1:N
            get_from_storage!(χ[k], X[k], n)
        end
        ϵₙ⁽ⁱ⁺¹⁾ = [ϵ⁽ⁱ⁾[l][n] for l ∈ 1:L]  # ϵₙ⁽ⁱ⁺¹⁾ ≈ ϵₙ⁽ⁱ⁾ for non-linear controls
        # TODO: we could add a self-consistent loop here for ϵₙ⁽ⁱ⁺¹⁾
        Δuₙ = zeros(L)  # Δu is Δϵ without (Sₗₙ/λₐ) factor
        for l = 1:L  # `l` is the index for the different controls
            for k = 1:N  # k is the index over the objectives
                ϕₖ = wrk.fw_propagators[k].state
                μₖₗ = wrk.control_derivs[k][l]
                if !isnothing(μₖₗ)
                    μₗₖₙ = _eval_mu(μₖₗ, wrk, ϵₙ⁽ⁱ⁺¹⁾, tlist, n)
                    Δuₙ[l] += Im(dot(χ[k], μₗₖₙ, ϕₖ))
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
            local ϕₖ = prop_step!(wrk.fw_propagators[k])
            write_to_storage!(Φ[k], n, ϕₖ)
        end
        # TODO: update sigma
    end  # time loop
end


function update_result!(wrk::KrotovWrk, i::Int64)
    res = wrk.result
    J_T_func = wrk.kwargs[:J_T]
    res.J_T_prev = res.J_T
    for k in eachindex(wrk.fw_propagators)
        res.states[k] = wrk.fw_propagators[k].state
    end
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
    for l in eachindex(ϵ_opt)
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
    flush(stdout)
end
