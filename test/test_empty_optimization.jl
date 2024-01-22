using Test
using StableRNGs
using QuantumControl: hamiltonian, optimize, ControlProblem, Trajectory
using QuantumControl.Controls: get_controls
using QuantumControlTestUtils.RandomObjects: random_matrix, random_state_vector

@testset "empty optimization" begin

    # Test that trying to run an optimization without any controls produces a
    # meaningful error message

    rng = StableRNG(2264511904)

    N = 10
    H = random_matrix(N; rng)
    trajectories = [
        Trajectory(
            random_state_vector(N; rng),
            H;
            target_state=random_state_vector(N; rng)
        )
    ]

    @test length(get_controls(trajectories)) == 0

    tlist = collect(range(0; length=1001, step=1.0))

    problem = ControlProblem(trajectories, tlist; pulse_options=Dict())

    msg = "no controls in trajectories: cannot optimize"
    @test_throws ErrorException(msg) optimize(problem; method=:krotov)

end
