using QuantumControlBase
using QuantumPropagators: initpropwrk, propstep!, init_storage, write_to_storage!, get_from_storage!
using LinearAlgebra
using Dates
import Base.Threads.@threads

"""Result object returned by [`optimize_pulses`](@ref)."""
mutable struct KrotovResult
    tlist :: Vector{Float64}
    iters :: Vector{Int64}
    iter_seconds :: Vector{Float64}
    tau_vals :: Vector{ComplexF64}
    guess_controls
    optimized_controls
    all_pulses
    states
    start_local_time
    end_local_time
    message

    function KrotovResult(problem)
        tlist = problem.tlist
        controls = getcontrols(problem.objectives)
        iters = Vector{Int64}()
        iter_seconds = Vector{Float64}()
        tau_vals = Vector{ComplexF64}()
        guess_controls = [
            discretize(control, tlist) for control in controls
        ]
        optimized_controls = [copy(guess) for guess in guess_controls]
        all_pulses = Vector{Any}()
        states = [similar(obj.initial_state) for obj in problem.objectives]
        start_local_time = now()
        end_local_time = now()
        message = "in progress"
        new(tlist, iters, iter_seconds, tau_vals, guess_controls,
            optimized_controls, all_pulses, states, start_local_time,
            end_local_time, message)
    end
end

Base.show(io::IO, r::KrotovResult) = print(io, "KrotovResult<$(r.message)>")
Base.show(io::IO, ::MIME"text/plain", r::KrotovResult) = print(io, """
Krotov Optimization Result
--------------------------
- Started at $(r.start_local_time)
- Number of objectives: $(length(r.states))
- Number of iterations: $(length(r.iters) - 1)
- Reason for termination: $(r.message)
- Ended at $(r.end_local_time) ($(r.end_local_time - r.start_local_time))""")


function update_result!(r::KrotovResult, i::Int64)
    push!(r.iters, i)
    prev_time = r.end_local_time
    r.end_local_time = now()
    iter_seconds = Dates.toms(r.end_local_time - prev_time) / 1000.0
    push!(r.iter_seconds, iter_seconds)
    # TODO
end

function finalize_result!(r::KrotovResult)
    r.end_local_time = now()
    r.message = "finished"  # TODO
    # TODO
end




# Krotov workspace (for internal use)
struct KrotovWrk
    # TODO: specify types more strictly

    # a copy of the objectives
    objectives :: Vector{QuantumControlBase.Objective}

    # the adjoint objectives, containing the adjoint generators for the
    # backward propagation
    adjoint_objectives :: Vector{QuantumControlBase.Objective}

    # The kwargs from the control problem
    kwargs :: AbstractDict

    # Tuple of the original controls (probably functions)
    controls :: Tuple

    # controls discretized on intervals of tlist
    guess_pulses :: Vector{Any}

    # The optimized pulses
    optimized_pulses :: Vector{Any}

    # map of controls to options
    pulse_options :: AbstractDict

    # Result object

    result :: KrotovResult

    #################################
    # scratch objects, per objective:

    # forward-propagated states
    # TODO: use result.states instead of ϕ
    ϕ :: Vector{Any}

    # backward-propagated states
    # TODO: rename to bw_states (maybe best not have unicode in field names)
    χ :: Vector{Any}

    # dynamical generator at a particular point in time
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
        kwargs = problem.kwargs
        guess_pulses = [
            discretize_on_midpoints(control, tlist) for control in controls
        ]
        optimized_pulses = [copy(pulse) for pulse in guess_pulses]
        pulse_options = problem.pulse_options
        result = KrotovResult(problem)
        ϕ = [similar(obj.initial_state) for obj in objectives]
        χ = [similar(obj.initial_state) for obj in objectives]
        zero_vals = IdDict(control => zero(guess_pulses[i][1]) for (i, control) in enumerate(controls))
        G = [setcontrolvals(obj.generator, zero_vals) for obj in objectives]
        vals_dict = [copy(zero_vals) for _ in objectives]
        # TODO: forward_storage only if g_b != 0
        fw_storage = [init_storage(obj.initial_state, tlist) for obj in objectives]
        # TODO: second forward storage only if second order
        fw_storage2 = [init_storage(obj.initial_state, tlist) for obj in objectives]
        bw_storage = [init_storage(obj.initial_state, tlist) for obj in objectives]
        prop_wrk = [initpropwrk(ϕ[i], tlist, G[i]) for i in 1:length(objectives)]
        # TODO: need separate prop_wrk for bw-prop?
        new(objectives, adjoint_objectives, kwargs, controls,
            guess_pulses, optimized_pulses, pulse_options, result, ϕ, χ, G,
            vals_dict, fw_storage, fw_storage2, bw_storage, prop_wrk)
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
keyword argument `chi!`.

# Optional problem keyword arguments

The following `problem` keyword arguments are supported (with default values):

* `sigma=nothing`: Function that calculate the second-order contribution. If
   not given, the first-order Krotov method is used.
* `iter_start=0`: the initial iteration number
* `iter_stop=5000`: the maximum iteration number

"""
function optimize_pulses(problem)
    #=chi! = problem.kwargs[:chi!]=#
    sigma = get(problem.kwargs, :sigma, nothing)
    iter_start = get(problem.kwargs, :iter_start, 0)
    iter_stop = get(problem.kwargs, :iter_stop, 5000)
    skip_initial_forward_propagation = get(
        problem.kwargs, :skip_initial_forward_propagation, false
    )

    wrk = KrotovWrk(problem)
    # TODO: if continuing previous optimization (`continue_from` argument),
    # ammend wrk from existing Result

    i = iter_start

    if skip_initial_forward_propagation
        @info "Skipping initial forward propagation"
    else
        @info "Initial Forward Propagation"
        @threads for (k, obj) in collect(enumerate(wrk.objectives))
            krotov_initial_fw_prop!(wrk.ϕ[k], obj.initial_state, k, wrk)
        end
    end

    # TODO: if sigma, fw_storage0 = fw_storage
    update_result!(wrk.result, 0)
    # TODO: info_hook(wrk)

    converged = (iter_stop <= i)
    # TODO: converged = converged || check_convergence(wrk)

    ϵ⁰ = wrk.guess_pulses
    ϵ¹ = wrk.optimized_pulses

    while !converged
        i = i + 1
        @info "Krotov iteration $i"
        krotov_iteration(wrk, ϵ⁰, ϵ¹)
        # TODO: info_hook(wrk, i)
        update_result!(wrk.result, i)
        converged = converged || (iter_stop <= i)
        # TODO: converged = converged || check_convergence(wrk)
        ϵ⁰, ϵ¹ = ϵ¹, ϵ⁰
    end

    finalize_result!(wrk.result)
    return wrk.result

end


function krotov_initial_fw_prop!(ϕₖ, ϕₖⁱⁿ, k, wrk)
    Φ₀ = wrk.fw_storage[k]
    copyto!(ϕₖ,  ϕₖⁱⁿ)
    (Φ₀ ≠ nothing) && write_to_storage!(Φ₀, 1,  ϕₖⁱⁿ)
    N_T = length(wrk.result.tlist) - 1
    ϵ⁰ = wrk.guess_pulses
    for n = 1:N_T
        G, dt = _fw_gen(ϵ⁰, k, n, wrk)
        propstep!(ϕₖ, G, dt, wrk.prop_wrk[k])
        (Φ₀ ≠ nothing) && write_to_storage!(Φ₀, n+1, ϕₖ)
    end
end


function krotov_iteration(wrk, ϵ⁰, ϵ¹)

    ϕ = wrk.ϕ  # assumed to contain the results of forward propagation
    χ = wrk.χ
    chi! = (χ, ϕ) -> copyto!(χ, ϕ)  # TODO: get function from arguments
    N_T = length(wrk.result.tlist) - 1
    N = length(wrk.objectives)
    L = length(wrk.controls)
    X = wrk.bw_storage
    Φ = wrk.fw_storage  # TODO: pass in Φ₁, Φ₀ as parameters

    # TODO: re-initialize wrk.prop_wrk?

    # backward propagation
    chi!(χ, ϕ)
    @threads for k = 1:N
        write_to_storage!(X[k], N_T+1, χ[k])
        for n = N_T:-1:1
            local (G, dt) = _bw_gen(ϵ⁰, k, n, wrk)
            propstep!(χ[k], G, dt, wrk.prop_wrk[k])
            write_to_storage!(X[k], n, χ[k])
        end
    end

    # pulse update and forward propagation

    @threads for k = 1:N
        copyto!(ϕ[k], wrk.objectives[k].initial_state)
    end

    for n = 1:N_T
        χ = wrk.χ
        for k = 1:N
            get_from_storage!(χ[k], X[k], n)
        end
        Δuₙ = zeros(L)
        for l = 1:L
            Sₗₙ_λₐₗ = 1.0 # TODO
            for k = 1:N
                μ_lkn = 1.0 # TODO
                Δuₙ[l] += Sₗₙ_λₐₗ * imag(dot(χ[k], μ_lkn, ϕ[k]))
            end
        end
        # TODO: second order update
        Δϵₙ = Δuₙ  # no parameterization (TODO)
        for l = 1:L
            ϵ¹[l][n] = ϵ⁰[l][n] + Δϵₙ[l]
        end
        @threads for k = 1:N
            local (G, dt) = _fw_gen(ϵ¹, k, n, wrk)
            propstep!(ϕ[k], G, dt, wrk.prop_wrk[k])
            write_to_storage!(Φ[k], n, ϕ[k])
        end
        # TODO: update sigma
    end
end


# The dynamical generator for the forward propagation
function _fw_gen(ϵ, k, n, wrk) # TODO
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
function _bw_gen(ϵ, k, n, wrk) # TODO: check if this is correct
    vals_dict = wrk.vals_dict[k]
    t = wrk.result.tlist
    for (l, control) in enumerate(wrk.controls)
        vals_dict[control] = ϵ[l][n]
    end
    dt = t[n+1] - t[n]
    setcontrolvals!(wrk.G[k], wrk.adjoint_objectives[k].generator, vals_dict)
    return wrk.G[k], -dt
end
