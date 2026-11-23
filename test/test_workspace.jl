using Test
using StableRNGs
using Krotov: KrotovWrk
using QuantumControl.DummyOptimization: dummy_control_problem
using QuantumControl.Functionals: J_T_ss
import QuantumPropagators
using QuantumPropagators: AbstractPropagator


# A propagation method that is not piecewise. The `method` is passed to
# `init_prop` as an instance of a type that we own, to avoid type piracy.
struct NonPiecewiseMethod end

struct NonPiecewisePropagator <: AbstractPropagator end

function QuantumPropagators.init_prop(
    state,
    generator,
    tlist,
    method::NonPiecewiseMethod;
    kwargs...
)
    return NonPiecewisePropagator()
end


function _problem(; kwargs...)
    rng = StableRNG(1244568944)
    return dummy_control_problem(;
        N = 2,
        density = 1.0,
        complex_operators = false,
        rng,
        J_T = J_T_ss,
        kwargs...
    )
end


function _test_rejects_non_piecewise(problem)
    @test_logs (:error, r"Cannot initialize propagation") match_mode = :any begin
        @test_throws "Cannot initialize propagator as piecewise" KrotovWrk(problem)
    end
end


@testset "piecewise propagators" begin

    wrk = KrotovWrk(_problem(prop_method = :cheby))
    @test wrk.fw_prop_kwargs[1][:piecewise] == true
    @test wrk.bw_prop_kwargs[1][:piecewise] == true
    @test wrk.bw_prop_kwargs[1][:backward] == true

    # forward propagation
    _test_rejects_non_piecewise(
        _problem(prop_method = :cheby, fw_prop_method = NonPiecewiseMethod())
    )
    # backward propagation
    _test_rejects_non_piecewise(
        _problem(prop_method = :cheby, bw_prop_method = NonPiecewiseMethod())
    )
    # The requirement cannot be disabled
    _test_rejects_non_piecewise(
        _problem(
            prop_method = NonPiecewiseMethod(),
            prop_piecewise = false,
            fw_prop_piecewise = false,
            bw_prop_piecewise = false,
        )
    )

    # The backward propagation cannot be disabled
    wrk = KrotovWrk(_problem(prop_method = :cheby, bw_prop_backward = false))
    @test wrk.bw_prop_kwargs[1][:backward] == true

end
