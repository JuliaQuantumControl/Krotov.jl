using Test
using TestItems

@testitem "iter_start_stop" begin
    # Test that setting iter_start and iter_stop in fact restricts the
    # optimization to those numbers

    using QuantumControl: optimize
    using StableRNGs
    using Krotov
    using QuantumControlTestUtils.DummyOptimization: dummy_control_problem
    using QuantumControl.Functionals: J_T_ss
    using IOCapture

    rng = StableRNG(1244568944)
    problem = dummy_control_problem(;
        iter_start=10,
        N=2,
        density=1.0,
        complex_operators=false,
        rng,
        J_T=J_T_ss,
        store_iter_info=["iter.", "J_T"]
    )
    captured = IOCapture.capture(passthrough=false) do
        optimize(problem; method=Krotov, iter_stop=12)
    end
    res = captured.value
    @test res.converged
    @test res.iter_start == 10
    @test res.iter_stop == 12
    iters = [values[1] for values in res.records]
    @test iters == [0, 11, 12]

end


@testitem "callback" begin

    using QuantumControl: optimize
    using LinearAlgebra: norm
    using LinearAlgebra.BLAS: scal!
    using Krotov
    using StableRNGs
    using QuantumControlTestUtils.DummyOptimization: dummy_control_problem
    using QuantumControl.Functionals: J_T_ss
    using IOCapture

    rng = StableRNG(1244568944)
    PASSTHROUGH = false


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
        N=2,
        density=1.0,
        complex_operators=false,
        rng,
        J_T=J_T_ss,
        callback=callback1,
    )

    captured = IOCapture.capture(passthrough=PASSTHROUGH) do
        optimize(problem; method=Krotov, iter_stop=1)
    end
    @test contains(
        captured.output,
        "This is callback 1 for iter 0\n iter.        J_T   ∫gₐ(t)dt          J       ΔJ_T         ΔJ    secs"
    )
    @test contains(captured.output, "This is callback 1 for iter 1\n     1")

    # passing `callback` to `optimize` overwrites `callback` in `problem`
    captured = IOCapture.capture(passthrough=PASSTHROUGH) do
        optimize(problem; method=Krotov, iter_stop=1, callback=callback2)
    end
    @test !contains(captured.output, "This is callback 1 for iter 0")
    @test !contains(captured.output, "This is callback 1 for iter 1")
    @test contains(captured.output, "This is callback 2 for iter 0")
    @test contains(captured.output, "This is callback 2 for iter 1")

    captured = IOCapture.capture(passthrough=PASSTHROUGH) do
        optimize(
            problem;
            method=Krotov,
            iter_stop=1,
            callback=(callback1, callback2),
            print_iters=false
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

    captured = IOCapture.capture(passthrough=PASSTHROUGH) do
        optimize(
            problem;
            method=Krotov,
            iter_stop=1,
            callback=(callback1, callback2),
            store_iter_info=["J_T"]
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
    captured = IOCapture.capture(passthrough=PASSTHROUGH) do
        optimize(
            problem;
            method=Krotov,
            iter_stop=3,
            callback=reduce_pulse,
            store_iter_info=["iter.", "J_T"]
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
