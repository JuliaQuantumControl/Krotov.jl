using Test
using QuantumControl: optimize
using StableRNGs
using LinearAlgebra: norm
using LinearAlgebra.BLAS: scal!
using Krotov
using QuantumControl.DummyOptimization: dummy_control_problem
using QuantumControl.Functionals: J_T_ss, make_chi
using QuantumPropagators.Controls: discretize
using IOCapture

PASSTHROUGH = false

@testset "iter_start_stop" begin
    # Test that setting iter_start and iter_stop in fact restricts the
    # optimization to those numbers
    rng = StableRNG(1244568944)
    problem = dummy_control_problem(;
        iter_start = 10,
        N = 2,
        density = 1.0,
        complex_operators = false,
        rng,
        J_T = J_T_ss,
        store_iter_info = ["iter.", "J_T"]
    )
    captured = IOCapture.capture(passthrough = PASSTHROUGH) do
        optimize(problem; method = Krotov, iter_stop = 12)
    end
    res = captured.value
    @test res.converged
    @test res.iter_start == 10
    @test res.iter_stop == 12
    iters = [values[1] for values in res.records]
    @test iters == [0, 11, 12]
end


@testset "callback" begin

    rng = StableRNG(1244568944)

    function callback1(_, iter, args...)
        println("This is callback 1 for iter $iter")
    end

    function callback2(_, iter, args...)
        println("This is callback 2 for iter $iter")
        return ("cb2", iter)
    end

    function reduce_pulse(wrk, iter, ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾)
        r0 = norm(ϵ⁽ⁱ⁾[1])
        r1 = norm(ϵ⁽ⁱ⁺¹⁾[1])
        scal!(0.8, ϵ⁽ⁱ⁺¹⁾[1])
        r2 = norm(ϵ⁽ⁱ⁺¹⁾[1])
        return (r0, r1, r2)
    end

    problem = dummy_control_problem(;
        N = 2,
        density = 1.0,
        complex_operators = false,
        rng,
        J_T = J_T_ss,
        callback = callback1,
    )

    captured = IOCapture.capture(passthrough = PASSTHROUGH) do
        optimize(problem; method = Krotov, iter_stop = 1)
    end
    @test contains(
        captured.output,
        "This is callback 1 for iter 0\n iter.        J_T   ∫gₐ(t)dt          J       ΔJ_T         ΔJ    secs"
    )
    @test contains(captured.output, "This is callback 1 for iter 1\n     1")

    # passing `callback` to `optimize` overwrites `callback` in `problem`
    captured = IOCapture.capture(passthrough = PASSTHROUGH) do
        optimize(problem; method = Krotov, iter_stop = 1, callback = callback2)
    end
    @test !contains(captured.output, "This is callback 1 for iter 0")
    @test !contains(captured.output, "This is callback 1 for iter 1")
    @test contains(captured.output, "This is callback 2 for iter 0")
    @test contains(captured.output, "This is callback 2 for iter 1")

    captured = IOCapture.capture(passthrough = PASSTHROUGH) do
        optimize(
            problem;
            method = Krotov,
            iter_stop = 1,
            callback = (callback1, callback2),
            print_iters = false
        )
    end
    @test captured.value.converged
    @test contains(
        captured.output,
        """
        This is callback 1 for iter 0
        This is callback 2 for iter 0
        This is callback 1 for iter 1
        This is callback 2 for iter 1
        """
    )
    @test captured.value.records == [("cb2", 0), ("cb2", 1)]

    captured = IOCapture.capture(passthrough = PASSTHROUGH) do
        optimize(
            problem;
            method = Krotov,
            iter_stop = 1,
            callback = (callback1, callback2),
            store_iter_info = ["J_T"]
        )
    end
    @test captured.value.converged
    @test length(captured.value.records) == 2
    @test length(captured.value.records[1]) == 3
    @test captured.value.records[1][1] == "cb2"
    @test captured.value.records[1][2] == 0
    @test captured.value.records[1][3] isa Float64

    # we should also be able to modify the updated pulses in the callback and
    # have that take effect.
    captured = IOCapture.capture(passthrough = PASSTHROUGH) do
        optimize(
            problem;
            method = Krotov,
            iter_stop = 3,
            callback = reduce_pulse,
            store_iter_info = ["iter.", "J_T"]
        )
    end
    @test captured.value.converged
    for i = 2:length(captured.value.records)
        record = captured.value.records[i]
        (nrm_guess, nrm_upd, nrm_upd_scaled, iter, J_T) = record
        nrm_upd_scaled_prev = captured.value.records[i-1][3]
        @test nrm_upd_scaled ≈ 0.8 * nrm_upd
        if i >= 3
            @test nrm_guess ≈ nrm_upd_scaled_prev
        end
    end

end


@testset "callback = nothing" begin
    # An explicit `callback = nothing` in the problem is equivalent to no
    # callback. We call the Krotov backend directly, since
    # `QuantumControl.optimize` may drop `callback = nothing`.
    rng = StableRNG(1244568944)
    problem = dummy_control_problem(;
        N = 2,
        density = 1.0,
        complex_operators = false,
        rng,
        J_T = J_T_ss,
        iter_stop = 1,
        callback = nothing,
    )
    res = optimize(problem, Val(:Krotov))
    @test res.converged
    @test res.iter == 1
    @test res.message == "Reached maximum number of iterations"
end


@testset "iter_stop reached before first iteration" begin
    # If `iter_start` (or the iteration of a `continue_from` result) already
    # reaches `iter_stop`, no iteration should be performed
    rng = StableRNG(1244568944)
    problem = dummy_control_problem(;
        N = 2,
        density = 1.0,
        complex_operators = false,
        rng,
        J_T = J_T_ss,
        print_iters = false,
    )

    res = optimize(problem; method = Krotov, iter_start = 5, iter_stop = 5)
    @test res.converged
    @test res.iter == 5
    @test res.message == "Reached maximum number of iterations"

    res1 = optimize(problem; method = Krotov, iter_stop = 3)
    @test res1.iter == 3
    controls1 = copy.(res1.optimized_controls)
    res2 = optimize(problem; method = Krotov, iter_stop = 3, continue_from = res1)
    @test res2 === res1  # updated in place
    @test res2.converged
    @test res2.iter == 3
    @test res2.message == "Reached maximum number of iterations"
    @test res2.optimized_controls ≈ controls1

    res3 = optimize(problem; method = Krotov, iter_stop = 5, continue_from = res1)
    @test res3.iter == 5
end


@testset "result describes last completed iteration" begin

    rng = StableRNG(1244568944)
    problem = dummy_control_problem(;
        N = 2,
        n_trajectories = 2,
        density = 1.0,
        complex_operators = false,
        rng,
        J_T = J_T_ss,
        print_iters = false,
    )
    N_T = length(problem.tlist) - 1  # number of time steps per propagation

    # Check that `res` is the result of an uninterrupted optimization to
    # iteration `iter`, and that it is consistent with its `states`
    function test_result_is_iteration(res, iter)
        @test res.iter == iter
        ref = optimize(problem; method = Krotov, iter_stop = iter)
        @test isapprox(res.J_T, J_T_ss(res.states, problem.trajectories); atol = 1e-12)
        @test isapprox(res.J_T, ref.J_T; atol = 1e-12)
        @test all(isapprox.(res.states, ref.states; atol = 1e-12))
        @test isapprox(res.tau_vals, ref.tau_vals; atol = 1e-12)
        @test all(isapprox.(res.optimized_controls, ref.optimized_controls; atol = 1e-12))
    end

    # A copy of the result during the optimization (as `@optimize_or_load`
    # writes when the process exits unexpectedly) can be continued
    snapshots = Dict{Int,Any}()
    function take_snapshot(wrk, iter, args...)
        snapshots[iter] = deepcopy(wrk.result)
        return nothing
    end
    res3 = optimize(problem; method = Krotov, iter_stop = 3, callback = take_snapshot)
    @test !(snapshots[2].optimized_controls ≈ snapshots[2].guess_controls)
    for iter = 1:3
        test_result_is_iteration(snapshots[iter], iter)
    end
    res2 = snapshots[2]
    states2 = res2.states
    state2 = states2[1]
    state2_copy = copy(state2)
    res2_continued = optimize(problem; method = Krotov, iter_stop = 3, continue_from = res2)
    @test res2_continued === res2  # updated in place
    test_result_is_iteration(res2_continued, 3)
    # The arrays of the previous iteration are not overwritten
    @test state2 !== res2.states[1]
    @test state2 == state2_copy

    # Interrupt the optimization in the n'th step of all propagations, which
    # must be a step of a backward propagation if `backward`, and of a forward
    # propagation otherwise. The initial forward propagation and the backward
    # propagation propagate the trajectories one after the other, the forward
    # propagation in each iteration propagates the trajectories together, step
    # by step.
    function interrupt_at(n; backward)
        n_steps = Ref(0)
        function interrupt(propagator, observables)
            n_steps[] += 1
            if n_steps[] == n
                @test propagator.backward == backward
                throw(InterruptException())
            end
            return nothing
        end
        return optimize(problem; method = Krotov, iter_stop = 5, prop_callback = interrupt)
    end
    # Before the guess is complete, the exception is not caught
    @test_throws InterruptException interrupt_at(N_T + N_T ÷ 2; backward = false)
    # backward propagation of iteration 2: first and second trajectory
    for n in (6 * N_T + N_T ÷ 2, 7 * N_T + N_T ÷ 2)
        res = interrupt_at(n; backward = true)
        @test startswith(res.message, "Exception: InterruptException")
        test_result_is_iteration(res, 1)
    end
    # forward propagation of iteration 2: first and second trajectory
    for n in (9 * N_T + 1, 9 * N_T + 2)
        res = interrupt_at(n; backward = false)
        @test startswith(res.message, "Exception: InterruptException")
        test_result_is_iteration(res, 1)
    end

    # Exception in `J_T` (in its third call, for iteration 2)
    n_calls = Ref(0)
    function raising_J_T(Ψ, trajectories)
        n_calls[] += 1
        (n_calls[] == 3) && throw(InterruptException())
        return J_T_ss(Ψ, trajectories)
    end
    res = optimize(
        problem;
        method = Krotov,
        iter_stop = 5,
        J_T = raising_J_T,
        chi = make_chi(J_T_ss, problem.trajectories)
    )
    @test startswith(res.message, "Exception: InterruptException")
    test_result_is_iteration(res, 1)

    # Exceptions in the `callback` or in `check_convergence`
    function raise_for_guess(wrk, iter, args...)
        (iter == 0) && throw(InterruptException())
        return nothing
    end
    @test_throws InterruptException begin
        optimize(problem; method = Krotov, iter_stop = 5, callback = raise_for_guess)
    end
    function raise_in_iter_2(wrk, iter, args...)
        (iter == 2) && throw(InterruptException())
        return nothing
    end
    res = optimize(problem; method = Krotov, iter_stop = 5, callback = raise_in_iter_2)
    @test startswith(res.message, "Exception: InterruptException")
    test_result_is_iteration(res, 2)
    function raise_at_2(result)
        (result.iter == 2) && throw(InterruptException())
        return result
    end
    res = optimize(problem; method = Krotov, iter_stop = 5, check_convergence = raise_at_2)
    @test startswith(res.message, "Exception: InterruptException")
    test_result_is_iteration(res, 2)

    # Pulses modified by a callback are reflected in `optimized_controls`,
    # which `check_convergence` receives after the callback
    modified_pulses = Vector{Float64}[]
    function reduce_pulse(wrk, iter, ϵ⁽ⁱ⁺¹⁾, ϵ⁽ⁱ⁾)
        if iter > 0
            scal!(0.8, ϵ⁽ⁱ⁺¹⁾[1])
            push!(modified_pulses, copy(ϵ⁽ⁱ⁺¹⁾[1]))
        end
        return nothing
    end
    seen_controls = Vector{Float64}[]
    function record_controls(result)
        (result.iter > 0) && push!(seen_controls, copy(result.optimized_controls[1]))
        return result
    end
    res = optimize(
        problem;
        method = Krotov,
        iter_stop = 2,
        callback = reduce_pulse,
        check_convergence = record_controls
    )
    @test length(seen_controls) == length(modified_pulses) == 2
    for (controls, pulse) in zip(seen_controls, modified_pulses)
        @test isapprox(controls, discretize(pulse, problem.tlist); atol = 1e-12)
    end
    @test isapprox(
        res.optimized_controls[1],
        discretize(modified_pulses[end], problem.tlist);
        atol = 1e-12
    )

end
