using Krotov:
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

function plot_positive_parametrization_comparison()

    u_vals = collect(range(-3, 3, length=101))
    ϵ_vals = collect(range(0, 1, length=101))
    ϵ_max = 1.0

    pnl1 = plot(
        u_vals,
        abs.(u_vals);
        linestyle=:dash,
        color="black",
        label="",
        xlabel="u",
        ylabel="ϵ",
        legend=false
    )
    plot!(pnl1, u_vals, TanhSqParametrization(ϵ_max).epsilon_of_u.(u_vals), label="TanhSq")
    plot!(
        pnl1,
        u_vals,
        LogisticSqParametrization(ϵ_max).epsilon_of_u.(u_vals),
        label="LogisticSq(k=1)"
    )
    plot!(
        pnl1,
        u_vals,
        LogisticSqParametrization(ϵ_max, k=4.0).epsilon_of_u.(u_vals),
        label="LogisticSq(k=4)"
    )
    plot!(pnl1, u_vals, SquareParametrization().epsilon_of_u.(u_vals), label="Square")
    ylims!(pnl1, (0, 1.2))

    pnl2 = plot(
        ϵ_vals,
        ϵ_vals;
        linestyle=:dash,
        color="black",
        label="",
        xlabel="ϵ",
        ylabel="u"
    )
    plot!(pnl2, ϵ_vals, TanhSqParametrization(ϵ_max).u_of_epsilon.(ϵ_vals), label="TanhSq")
    plot!(
        pnl2,
        ϵ_vals,
        LogisticSqParametrization(ϵ_max).u_of_epsilon.(ϵ_vals),
        label="LogisticSq(k=1)"
    )
    plot!(
        pnl2,
        ϵ_vals,
        LogisticSqParametrization(ϵ_max, k=4.0).u_of_epsilon.(ϵ_vals),
        label="LogisticSq(k=4)"
    )
    plot!(pnl2, ϵ_vals, SquareParametrization().u_of_epsilon.(ϵ_vals), label="Square")
    ylims!(pnl2, (0, 3))

    pnl3 = plot(
        u_vals,
        sign.(u_vals);
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
        TanhSqParametrization(ϵ_max).de_du_derivative.(u_vals),
        label="TanhSq"
    )
    plot!(
        pnl3,
        u_vals,
        LogisticSqParametrization(ϵ_max).de_du_derivative.(u_vals),
        label="LogisticSq(k=1)"
    )
    plot!(
        pnl3,
        u_vals,
        LogisticSqParametrization(ϵ_max, k=4.0).de_du_derivative.(u_vals),
        label="LogisticSq(k=4)"
    )
    plot!(pnl3, u_vals, SquareParametrization().de_du_derivative.(u_vals), label="Square")
    ylims!(pnl3, (-2, 2))

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
    gui(plot_positive_parametrization_comparison())
    readline()
end
