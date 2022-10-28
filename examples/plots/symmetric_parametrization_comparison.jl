using QuantumControl.PulseParametrizations:
    SquareParametrization,
    TanhParametrization,
    TanhSqParametrization,
    LogisticParametrization,
    LogisticSqParametrization

using Plots
Plots.default(
    linewidth               = 3,
    size                    = (550, 300),
    legend                  = :top,
    foreground_color_legend = nothing,
    background_color_legend = RGBA(1, 1, 1, 0.8),
)

function plot_symmetric_parametrization_comparison()

    u_vals = collect(range(-3, 3, length=101))
    ϵ_vals = collect(range(-1, 1, length=101))

    ϵ_min = -1.0
    ϵ_max = 1.0

    pnl1 = plot(
        u_vals,
        u_vals;
        linestyle=:dash,
        color="black",
        label="",
        xlabel="u",
        ylabel="ϵ",
        legend=false
    )
    plot!(
        pnl1,
        u_vals,
        TanhParametrization(ϵ_min, ϵ_max).a_of_epsilon.(u_vals),
        label="Tanh"
    )
    plot!(
        pnl1,
        u_vals,
        LogisticParametrization(ϵ_min, ϵ_max).a_of_epsilon.(u_vals),
        label="Logistic(k=1)"
    )
    plot!(
        pnl1,
        u_vals,
        LogisticParametrization(ϵ_min, ϵ_max, k=4).a_of_epsilon.(u_vals),
        label="Logistic(k=4)"
    )
    ylims!(pnl1, (-1.2, 1.2))

    pnl2 = plot(
        ϵ_vals,
        ϵ_vals;
        linestyle=:dash,
        color="black",
        label="",
        xlabel="ϵ",
        ylabel="u"
    )
    plot!(
        pnl2,
        ϵ_vals,
        TanhParametrization(ϵ_min, ϵ_max).epsilon_of_a.(ϵ_vals),
        label="Tanh"
    )
    plot!(
        pnl2,
        ϵ_vals,
        LogisticParametrization(ϵ_min, ϵ_max).epsilon_of_a.(ϵ_vals),
        label="Logistic(k=1)"
    )
    plot!(
        pnl2,
        ϵ_vals,
        LogisticParametrization(ϵ_min, ϵ_max, k=4).epsilon_of_a.(ϵ_vals),
        label="Logistic(k=4)"
    )
    ylims!(pnl2, (-3, 3))

    pnl3 = plot(
        u_vals,
        [1.0 for _ in u_vals];
        linestyle=:dash,
        color="black",
        label="",
        xlabel="u",
        ylabel="∂ϵ/∂u",
        legend=false
    )
    plot!(
        pnl3,
        u_vals,
        TanhParametrization(ϵ_min, ϵ_max).da_deps_derivative.(u_vals),
        label="Tanh"
    )
    plot!(
        pnl3,
        u_vals,
        LogisticParametrization(ϵ_min, ϵ_max).da_deps_derivative.(u_vals),
        label="Logistic(k=1)"
    )
    plot!(
        pnl3,
        u_vals,
        LogisticParametrization(ϵ_min, ϵ_max, k=4).da_deps_derivative.(u_vals),
        label="Logistic(k=4)"
    )
    ylims!(pnl3, (0, 2))

    plot(
        pnl1,
        pnl2,
        pnl3,
        layout=(1, 3),
        size=(1000, 300),
        left_margin=20Plots.px,
        bottom_margin=20Plots.px
    )

end

if abspath(PROGRAM_FILE) == @__FILE__
    gui(plot_symmetric_parametrization_comparison())
    readline()
end
