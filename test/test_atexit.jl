using Test
using StableRNGs
using QuantumControl: optimize, load_optimization, remove_atexit_save_optimization
using QuantumControl.DummyOptimization: dummy_control_problem
using QuantumControl.Functionals: J_T_ss
using Krotov


# The `atexit_filename` keyword argument (passed by `@optimize_or_load`)
# registers a hook that writes the result to the file when Julia exits. The
# hook must be removed when `optimize` returns or throws an exception, or it
# would overwrite the file at exit.
@testset "atexit hook" begin

    rng = StableRNG(1244568944)
    problem = dummy_control_problem(;
        N = 2,
        density = 1.0,
        complex_operators = false,
        rng,
        J_T = J_T_ss,
        print_iters = false,
    )
    atexit_filename = joinpath(mktempdir(), "result.jld2")

    # Compile `optimize` before taking a snapshot of `Base.atexit_hooks`, in
    # case anything registers a hook when it first runs
    optimize(problem; method = Krotov, iter_stop = 1, atexit_filename)
    hooks = copy(Base.atexit_hooks)
    is_new(hook) = !any(h -> h === hook, hooks)
    new_hooks() = filter(is_new, Base.atexit_hooks)

    # Check that there are no new hooks, and remove any leaked hooks so that a
    # failure does not affect later tests
    function test_no_new_hooks()
        leaked = new_hooks()
        @test isempty(leaked)
        foreach(remove_atexit_save_optimization, leaked)
    end

    function raise_in_iter_2(wrk, iter, args...)
        (iter == 2) && error("failure in iteration 2")
        return nothing
    end

    # The hook is registered during the optimization, and writes the current
    # result
    hooks_during_optimization = []
    function run_new_hooks(wrk, iter, args...)
        if iter == 1
            append!(hooks_during_optimization, new_hooks())
            foreach(hook -> hook(), hooks_during_optimization)
        end
        return nothing
    end
    res = optimize(
        problem;
        method = Krotov,
        iter_stop = 2,
        callback = run_new_hooks,
        atexit_filename
    )
    @test length(hooks_during_optimization) == 1
    saved = load_optimization(atexit_filename)
    @test saved.iter == 1
    @test saved.message == "Abort: ATEXIT"
    @test res.iter == 2
    test_no_new_hooks()

    # normal return
    res = optimize(problem; method = Krotov, iter_stop = 2, atexit_filename)
    @test res.iter == 2
    test_no_new_hooks()

    # exception that is caught
    res = optimize(
        problem;
        method = Krotov,
        iter_stop = 5,
        callback = raise_in_iter_2,
        atexit_filename
    )
    @test startswith(res.message, "Exception: failure in iteration 2")
    test_no_new_hooks()

    # exception that is rethrown
    @test_throws "failure in iteration 2" begin
        optimize(
            problem;
            method = Krotov,
            iter_stop = 5,
            callback = raise_in_iter_2,
            rethrow_exceptions = true,
            atexit_filename
        )
    end
    test_no_new_hooks()

    # another hook registered during the optimization stays registered
    other_hook() = nothing
    function register_other_hook(wrk, iter, args...)
        (iter == 1) && atexit(other_hook)
        return nothing
    end
    res = optimize(
        problem;
        method = Krotov,
        iter_stop = 2,
        callback = register_other_hook,
        atexit_filename
    )
    @test new_hooks() == [other_hook]
    filter!(h -> h !== other_hook, Base.atexit_hooks)
    test_no_new_hooks()

end


# When Julia exits, the hook writes the result only if the optimization is
# still running
@testset "atexit hook at exit" begin

    dir = mktempdir()
    script = joinpath(dir, "script.jl")
    write(
        script,
        raw"""
        using QuantumControl: optimize
        using QuantumControl.DummyOptimization: dummy_control_problem
        using QuantumControl.Functionals: J_T_ss
        using StableRNGs
        using Krotov

        dir = ARGS[1]
        problem = dummy_control_problem(;
            N = 2,
            density = 1.0,
            complex_operators = false,
            rng = StableRNG(1244568944),
            J_T = J_T_ss,
            print_iters = false,
        )

        # An optimization that ends with an exception must not write its file
        function raise_in_iter_2(wrk, iter, args...)
            (iter == 2) && error("failure in iteration 2")
            return nothing
        end
        try
            optimize(
                problem;
                method = Krotov,
                iter_stop = 5,
                callback = raise_in_iter_2,
                rethrow_exceptions = true,
                atexit_filename = joinpath(dir, "rethrown.jld2")
            )
        catch
        end

        # An optimization that is running when Julia exits writes its file
        function exit_in_iter_2(wrk, iter, args...)
            (iter == 2) && exit(3)
            return nothing
        end
        optimize(
            problem;
            method = Krotov,
            iter_stop = 5,
            callback = exit_in_iter_2,
            atexit_filename = joinpath(dir, "exited.jld2")
        )
        """
    )
    cmd =
        `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $script $dir`
    proc = run(pipeline(ignorestatus(cmd); stdout = devnull, stderr = devnull))
    @test proc.exitcode == 3
    @test !isfile(joinpath(dir, "rethrown.jld2"))
    @test isfile(joinpath(dir, "exited.jld2"))
    if isfile(joinpath(dir, "exited.jld2"))
        saved = load_optimization(joinpath(dir, "exited.jld2"))
        @test saved.iter == 2  # the `callback` runs after iteration 2
        @test saved.message == "Abort: ATEXIT"
    end

end
