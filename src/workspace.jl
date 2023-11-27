import QuantumControlBase
using QuantumControlBase: get_control_derivs
using QuantumControlBase.QuantumPropagators.Controls: get_controls, discretize_on_midpoints
using QuantumControlBase.QuantumPropagators: init_prop
using QuantumControlBase.QuantumPropagators.Storage: init_storage
using ConcreteStructs

# Krotov workspace (for internal use)
@concrete terse struct KrotovWrk

    # a copy of the objectives
    objectives

    # the adjoint objectives, containing the adjoint generators for the
    # backward propagation
    adjoint_objectives

    # The kwargs from the control problem
    kwargs

    # Tuple of the original controls (probably functions)
    controls

    # storage for controls discretized on intervals of tlist
    pulses0::Vector{Vector{Float64}}

    # second pulse storage: pulses0 and pulses1 alternate in storing the guess
    # pulses and optimized pulses in each iteration
    pulses1::Vector{Vector{Float64}}

    # values of ∫gₐ(t)dt for each pulse
    g_a_int::Vector{Float64}

    # update shapes S(t) for each pulse, discretized on intervals
    update_shapes::Vector{Vector{Float64}}

    lambda_vals::Vector{Float64}

    # map of controls to options
    pulse_options

    # Result object

    result

    #################################
    # scratch objects, per objective:

    control_derivs

    fw_storage # forward storage array (per objective)

    fw_storage2 # forward storage array (per objective)

    bw_storage # backward storage array (per objective)

    fw_propagators

    bw_propagators

    use_threads::Bool

end


function KrotovWrk(problem::QuantumControlBase.ControlProblem; verbose=false)
    use_threads = get(problem.kwargs, :use_threads, false)
    objectives = [obj for obj in problem.objectives]
    adjoint_objectives = [adjoint(obj) for obj in problem.objectives]
    controls = get_controls(objectives)
    if length(controls) == 0
        error("no controls in objectives: cannot optimize")
    end
    control_derivs = [get_control_derivs(obj.generator, controls) for obj in objectives]
    tlist = problem.tlist
    kwargs = Dict(problem.kwargs)  # creates a shallow copy; ok to modify
    default_update_shape = get(problem.kwargs, :update_shape, t -> 1.0)
    default_lambda_a = convert(Float64, get(problem.kwargs, :lambda_a, 1.0))
    default_pulse_options = IdDict(
        c => Dict(:lambda_a => default_lambda_a, :update_shape => default_update_shape)
        for c ∈ controls
    )
    if :pulse_options in keys(kwargs)
        if :update_shape in keys(kwargs)
            @warn("`update_shape` is ignored due to given `pulse_options`")
        end
        if :lambda_a in keys(kwargs)
            @warn("`lambda_a=$lambda_a` is ignored due to given `pulse_options`")
        end
    else
        if (:update_shape ∉ keys(kwargs)) && (:lambda_a ∉ keys(kwargs))
            @warn "Using default pulse_options: (:lambda_a => 1.0, :update_shape => (t -> 1.0))"
        end
    end
    pulse_options = get(kwargs, :pulse_options, default_pulse_options)
    for c ∈ controls
        if c ∉ keys(pulse_options)
            error("pulse_options must be defined for all controls")
        end
    end
    update_shapes = [
        discretize_on_midpoints(pulse_options[control][:update_shape], tlist) for
        control in controls
    ]
    lambda_vals =
        [convert(Float64, pulse_options[control][:lambda_a]) for control in controls]
    if haskey(kwargs, :continue_from)
        @info "Continuing previous optimization"
        result = kwargs[:continue_from]
        if !(result isa KrotovResult)
            # account for continuing from a different optimization method
            result = convert(KrotovResult, result)
        end
        result.iter_stop = get(kwargs, :iter_stop, 5000)
        result.converged = false
        result.start_local_time = now()
        result.message = "in progress"
        pulses0 = [
            discretize_on_midpoints(control, tlist) for control in result.optimized_controls
        ]
    else
        result = KrotovResult(problem)
        pulses0 = [discretize_on_midpoints(control, tlist) for control in controls]
    end
    pulses1 = [copy(pulse) for pulse in pulses0]
    g_a_int = zeros(length(pulses0))
    # TODO: forward_storage only if g_b != 0
    fw_storage = [init_storage(obj.initial_state, tlist) for obj in objectives]
    # TODO: second forward storage only if second order
    fw_storage2 = [init_storage(obj.initial_state, tlist) for obj in objectives]
    bw_storage = [init_storage(obj.initial_state, tlist) for obj in objectives]
    kwargs[:piecewise] = true  # only accept piecewise propagators
    fw_prop_method = [
        QuantumControlBase.get_objective_prop_method(
            obj,
            :fw_prop_method,
            :prop_method;
            kwargs...
        ) for obj in objectives
    ]
    bw_prop_method = [
        QuantumControlBase.get_objective_prop_method(
            obj,
            :bw_prop_method,
            :prop_method;
            kwargs...
        ) for obj in objectives
    ]

    fw_propagators = [
        begin
            verbose &&
                @info "Initializing fw-prop of objective $k with method $(fw_prop_method[k])"
            init_prop(
                obj.initial_state,
                obj.generator,
                tlist;
                method=fw_prop_method[k],
                kwargs...
            )
        end for (k, obj) in enumerate(objectives)
    ]
    bw_propagators = [
        begin
            verbose &&
                @info "Initializing bw-prop of objective $k with method $(bw_prop_method[k])"
            init_prop(
                obj.initial_state,
                obj.generator,
                tlist;
                method=bw_prop_method[k],
                backward=true,
                kwargs...
            )
        end for (k, obj) in enumerate(adjoint_objectives)
    ]
    return KrotovWrk(
        objectives,
        adjoint_objectives,
        kwargs,
        controls,
        pulses0,
        pulses1,
        g_a_int,
        update_shapes,
        lambda_vals,
        pulse_options,
        result,
        control_derivs,
        fw_storage,
        fw_storage2,
        bw_storage,
        fw_propagators,
        bw_propagators,
        use_threads
    )
end
