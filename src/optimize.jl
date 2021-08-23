using QuantumControlBase
using QuantumPropagators
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

    fw_storage :: Vector{Any}

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
        fw_storage = [init_storage(obj.initial_state, tlist) for obj in objectives]
        prop_wrk = [initpropwrk(state[i], tlist, G[i]) for i in 1:length(objectives)]
        new(objectives, adjoint_objectives, tlist, controls,
            guess_pulses, pulses, pulse_options, state, G, vals_dict,
            fw_storage, prop_wrk)
    end

end


function optimize_pulses(
    control_problem;
    #=chi_constructor,=#
    # TODO: use kwargs from ControlProblem
    sigma=nothing,
    iter_start=0,
    iter_stop=5000,
    skip_initial_forward_propagation=false,
)
    wrk = KrotovWrk(control_problem)

    if skip_initial_forward_propagation
        @info "Skipping initial forward propagation"
    else
        @info "Initial Forward Propagation"
        @threads for (i_obj, obj) in collect(enumerate(wrk.objectives))
            _krotov_initial_fw_prop!(
                wrk.state[i_obj], obj.initial_state, wrk.tlist,
                wrk.controls, wrk.vals_dict[i_obj], wrk.guess_pulses,
                wrk.G[i_obj], wrk.objectives[i_obj].generator,
                wrk.prop_wrk[i_obj], wrk.fw_storage[i_obj]
            )
        end
    end

    return wrk # XXX

end


function _krotov_initial_fw_prop!(Ψ, Ψ₀, tlist, controls, vals_dict,
                                  guess_pulses, G, generator, prop_wrk, storage)
    copyto!(Ψ, Ψ₀)
    write_to_storage!(storage, 1, Ψ₀)
    nsteps = length(tlist) - 1
    for i = 1:nsteps
        for (i_ctrl, control) in enumerate(controls)
            vals_dict[control] = guess_pulses[i_ctrl][i]
        end
        dt = tlist[i+1] - tlist[i]
        setcontrolvals!(G, generator, vals_dict)
        QuantumPropagators.propstep!(Ψ, G, dt, prop_wrk)
        write_to_storage!(storage, i+1, Ψ)
    end
end
