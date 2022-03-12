using QuantumControlBase.QuantumPropagators: propstep!, write_to_storage!, get_from_storage!
using QuantumControlBase: evalcontrols!
using QuantumControlBase.Functionals: make_chi
using QuantumControlBase.ConditionalThreads: @threadsif
using LinearAlgebra
using Printf


"""Optimize a control problem using Krotov's method.

```julia
result = optimize_krotov(problem)
```

optimizes the given
control [`problem`](@ref QuantumControlBase.ControlProblem),
returning a [`KrotovResult`](@ref).

!!! note

    It is recommended to call [`optimize`](@ref QuantumControlBase.optimize)
    with `method=:krotov` instead of calling `optimize_krotov` directly.

Keyword arguments that control the optimization are taken from the keyword
arguments used in the instantiation of `problem`.

# Required problem keyword arguments

* `J_T`: A function `J_T(ϕ, objectives)` that evaluates the final time
  functional from a list `ϕ` of forward-propagated states and
  `problem.objectives`.

# Optional problem keyword arguments

The following keyword arguments are supported (with default values):

* `chi`: A function `chi!(χ, ϕ, objectives)` what receives a list `ϕ`
  of the forward propagates state and must set ``-χₖ=∂J_T/∂⟨ϕₖ|``. If not
  given, it will be automatically determined from `J_T` via [`make_chi`](@ref)
* `force_zygote=false`: Whether to force the use of automatic differentiation
  when calling [`make_chi`](@ref).
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
function optimize_krotov(problem)
    sigma = get(problem.kwargs, :sigma, nothing)
    iter_start = get(problem.kwargs, :iter_start, 0)
    update_hook! = get(problem.kwargs, :update_hook, (args...) -> nothing)
    info_hook = get(problem.kwargs, :info_hook, print_table)
    check_convergence! = get(problem.kwargs, :check_convergence, res -> res)
    # note: the default `check_convergence!` is a no-op. We still always check
    # for "Reached maximum number of iterations" in `update_result!`
    skip_initial_forward_propagation =
        get(problem.kwargs, :skip_initial_forward_propagation, false)

    wrk = KrotovWrk(problem)

    ϵ⁽ⁱ⁾ = wrk.pulses0
    ϵ⁽ⁱ⁺¹⁾ = wrk.pulses1

    if skip_initial_forward_propagation
        @info "Skipping initial forward propagation"
    else
        @threadsif wrk.use_threads for (k, obj) in collect(enumerate(wrk.objectives))
            krotov_initial_fw_prop!(ϵ⁽ⁱ⁾, wrk.result.states[k], obj.initial_state, k, wrk)
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
    copyto!(ϕₖ, ϕₖⁱⁿ)
    (Φ₀ !== nothing) && write_to_storage!(Φ₀, 1, ϕₖⁱⁿ)
    N_T = length(wrk.result.tlist) - 1
    for n = 1:N_T
        G, dt = _fw_gen(ϵ⁽⁰⁾, k, n, wrk)
        propstep!(ϕₖ, G, dt, wrk.fw_prop_wrk[k])
        (Φ₀ !== nothing) && write_to_storage!(Φ₀, n + 1, ϕₖ)
    end
    # TODO: allow a custom propstep! routine
end


function krotov_iteration(wrk, ϵ⁽ⁱ⁾, ϵ⁽ⁱ⁺¹⁾)

    ϕ = wrk.result.states  # assumed to contain the results of forward propagation
    χ = wrk.bw_states
    J_T_func = wrk.kwargs[:J_T]
    force_zygote = get(wrk.kwargs, :get_zygote, false)
    chi! = get(wrk.kwargs, :chi, make_chi(J_T_func, wrk.objectives; force_zygote))
    N_T = length(wrk.result.tlist) - 1
    N = length(wrk.objectives)
    L = length(wrk.controls)
    X = wrk.bw_storage
    Φ = wrk.fw_storage  # TODO: pass in Φ₁, Φ₀ as parameters
    ∫gₐdt = wrk.g_a_int
    Im = imag

    # backward propagation
    chi!(χ, ϕ, wrk.objectives)
    @threadsif wrk.use_threads for k = 1:N
        write_to_storage!(X[k], N_T + 1, χ[k])
        for n = N_T:-1:1
            local (G, dt) = _bw_gen(ϵ⁽ⁱ⁾, k, n, wrk)
            propstep!(χ[k], G, dt, wrk.bw_prop_wrk[k])
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
        ϵₙ⁽ⁱ⁺¹⁾ = [ϵ⁽ⁱ⁾[l][n] for l ∈ 1:L]  # ϵₙ⁽ⁱ⁺¹⁾ ≈ ϵₙ⁽ⁱ⁾ for non-linear controls
        # TODO: we could add a self-consistent loop here for ϵₙ⁽ⁱ⁺¹⁾
        Δuₙ′ = zeros(L)  # for step size 1
        for l = 1:L  # `l` is the index for the different controls
            ∂ϵₗ╱∂u::Function = wrk.parametrization[l].de_du_derivative
            uₗ = wrk.parametrization[l].u_of_epsilon
            for k = 1:N  # k is the index over the objectives
                ∂Hₖ╱∂ϵₗ::Union{Function,Nothing} = wrk.control_derivs[k][l]
                if !isnothing(∂Hₖ╱∂ϵₗ)
                    μₗₖₙ = (∂Hₖ╱∂ϵₗ)(ϵₙ⁽ⁱ⁺¹⁾[l])
                    if wrk.is_parametrized[l]
                        ∂ϵₗ╱∂uₗ = (∂ϵₗ╱∂u)(uₗ(ϵₙ⁽ⁱ⁺¹⁾[l]))
                        Δuₙ′[l] += ∂ϵₗ╱∂uₗ * Im(dot(χ[k], μₗₖₙ, ϕ[k]))
                    else
                        Δuₙ′[l] += Im(dot(χ[k], μₗₖₙ, ϕ[k]))
                    end
                end
            end
        end
        # TODO: second order update
        for l = 1:L
            Sₗ = wrk.update_shapes[l]
            λₐ = wrk.lambda_vals[l]
            αₗ = (Sₗ[n] / λₐ)  # Krotov step size
            Δuₗₙ = αₗ * Δuₙ′[l]
            if wrk.is_parametrized[l]
                uₗ = wrk.parametrization[l].u_of_epsilon
                ϵₗ = wrk.parametrization[l].epsilon_of_u
                ϵ⁽ⁱ⁺¹⁾[l][n] = ϵₗ(uₗ(ϵ⁽ⁱ⁾[l][n]) + Δuₗₙ)
            else
                Δϵₗₙ = Δuₗₙ
                ϵ⁽ⁱ⁺¹⁾[l][n] = ϵ⁽ⁱ⁾[l][n] + Δϵₗₙ
            end
            (∫gₐdt)[l] += αₗ * abs(Δuₙ′[l])^2 * dt
        end
        # TODO: end of self-consistent loop
        @threadsif wrk.use_threads for k = 1:N
            local (G, dt) = _fw_gen(ϵ⁽ⁱ⁺¹⁾, k, n, wrk)
            propstep!(ϕ[k], G, dt, wrk.fw_prop_wrk[k])
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
    for l = 1:length(ϵ_opt)
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
