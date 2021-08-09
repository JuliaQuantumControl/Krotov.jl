"""Krotov Workspace."""
struct KrotovWrk
    # copy of objectives
    # adjoint adjectives
    # pulse mapping
    # Temporary Hamiltonians
    # Temporary States
    # Storage Arrays
end


function optimize_pulses(
    objectives,
    pulse_options,
    tlist;
    propagator,
    chi_constructor,
    sigma=nothing,
    iter_start=0,
    iter_stop=5000,
)
    wrk = KrotovWrk(objectives) # TODO

    # Initial forward propagation
    for (i_obj, obj) in enumerate(objectives)
        copyto!(wrk.state[i_obj], obj.initial_state)
        set_controls(wrk.H[i_obj], obj.H, pulse_vals)
        # TODO
    end

end
