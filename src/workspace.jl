import QuantumControlBase
using QuantumControlBase: getcontrols, getcontrolderivs, discretize_on_midpoints, evalcontrols
using QuantumPropagators: init_storage

# Krotov workspace (for internal use)
struct KrotovWrk{
        OT<:QuantumControlBase.AbstractControlObjective,
        AOT<:QuantumControlBase.AbstractControlObjective,
        KWT, CTRST<:Tuple, POT<:AbstractDict, STST, VDT, STORT, FWPRWT, BWPRWT,
        GT
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

    fw_prop_wrk :: Vector{FWPRWT}

    bw_prop_wrk :: Vector{BWPRWT}

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
        fw_prop_wrk = [
            QuantumControlBase.initobjpropwrk(obj, tlist, prop_method;
                                              initial_state=obj.initial_state,
                                              kwargs...)
            for obj in objectives
        ]
        bw_prop_wrk = [
            QuantumControlBase.initobjpropwrk(obj, tlist, prop_method;
                                              initial_state=obj.initial_state,
                                              kwargs...)
            for obj in objectives
        ]
        new{eltype(objectives),eltype(adjoint_objectives), typeof(kwargs),
            typeof(controls), typeof(pulse_options),
            typeof(objectives[1].initial_state), eltype(vals_dict),
            eltype(fw_storage), eltype(fw_prop_wrk), eltype(bw_prop_wrk),
            eltype(G)
        }(
            objectives, adjoint_objectives, kwargs, controls, pulses0, pulses1,
            g_a_int, update_shapes, lambda_vals, is_parametrized,
            parametrization, pulse_options, result, bw_states, G,
            control_derivs, vals_dict, fw_storage, fw_storage2,
            bw_storage, fw_prop_wrk, bw_prop_wrk, use_threads
        )
    end

end
