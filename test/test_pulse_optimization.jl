using Test
using QuantumControl
using LinearAlgebra
using StableRNGs
using QuantumControlTestUtils.DummyOptimization: dummy_control_problem
using QuantumControl.Controls: get_controls, discretize_on_midpoints
using QuantumControl.Functionals: J_T_re

@testset "pulse optimization" begin

    # Test the resolution of #28

    rng = StableRNG(1244561944)

    # The problem occurs when the controls are actually pulses (on the
    # midpoints of the time grid), so that the optimization does not have to
    # call `discretize_on_midpoints` internally
    problem = dummy_control_problem(; pulses_as_controls=true)
    nt = length(problem.tlist)
    guess_pulse = QuantumControl.Controls.get_controls(problem.trajectories)[1]
    @test length(guess_pulse) == nt - 1
    guess_pulse_copy = copy(QuantumControl.Controls.get_controls(problem.trajectories)[1])

    # Optimizing this should not modify the original generator in any way
    res = optimize(problem; method=:krotov, J_T=J_T_re, iter_stop=2)
    opt_control = res.optimized_controls[1]
    @test length(opt_control) == nt  # optimized_controls are always *on* tlist
    opt_pulse = discretize_on_midpoints(opt_control, problem.tlist)
    post_pulse = QuantumControl.Controls.get_controls(problem.trajectories)[1]

    # * The generator should still have the exact same objects as controls
    @test guess_pulse ≡ post_pulse
    # * These objects should not have been modified
    @test norm(guess_pulse_copy - guess_pulse) ≈ 0.0
    # * But the values of the optimized pulse should differ from the pulse in
    #   the generator
    @test norm(post_pulse - opt_pulse) > 0.1

    # The underlying reason for why this might fail is if
    # `discretize_on_midpoints` does not create a copy
    @test discretize_on_midpoints(guess_pulse, problem.tlist) ≢ guess_pulse

end
