"""Specification for a "time-local" pulse parametrization.

The parametrization is given as a collection of three functions:

* ``ϵ(u(t))``
* ``u(ϵ(t))``
* ``∂ϵ/∂u`` as a function of ``u(t)``.
"""
struct PulseParametrization
    epsilon_of_u :: Function
    u_of_epsilon :: Function
    de_du_derivative :: Function
end


"""Parametrization where ``ϵ(t) ≡ u(t)``."""
NoParametrization() = PulseParametrization(u->u, ϵ->ϵ, u->1.0)


"""Parametrization ϵ(t) = u²(t), enforcing pulse values ≥ 0."""
SquareParametrization() = PulseParametrization(
    u -> begin # epsilon_of_u
        ϵ = u^2
    end,
    ϵ -> begin # u_of_epsilon
        a = max(ϵ, 0.0)
        u = √a
    end,
    u -> begin # de_du_derivative
        ∂ϵ╱∂u = 2u
    end
)


"""Parametrization with a tanh function that enforces ϵ_min < ϵ(t) < ϵ_max.
"""
function TanhParametrization(ϵ_min, ϵ_max)

    Δ = ϵ_max - ϵ_min
    Σ = ϵ_max + ϵ_min
    εₚ = eps(1.0)  # 2⋅10⁻¹⁶ (machine precision)
    @assert ϵ_max > ϵ_min

    PulseParametrization(
        u -> begin # epsilon_of_u
            ϵ = tanh(u) * Δ / 2 + Σ / 2
        end,
        ϵ -> begin # u_of_epsilon
            a = clamp(2ϵ/Δ -  Σ/Δ, -1+εₚ, 1-εₚ)
            u = atanh(a)  # -18.4 < u < 18.4
        end,
        u -> begin # de_du_derivative
            ∂ϵ╱∂u = (Δ / 2) * sech(u)^2
        end
    )

end


"""Parametrization with a tanh² function that enforces ``0 ≤ ϵ(t) < ϵ_max``.
"""
function TanhSqParametrization(ϵ_max)

    εₚ = eps(1.0)  # 2⋅10⁻¹⁶ (machine precision)
    @assert ϵ_max > 0

    PulseParametrization(
        u -> begin # epsilon_of_u
            ϵ = ϵ_max * tanh(u)^2
        end,
        ϵ -> begin # u_of_epsilon
            a = clamp(ϵ / ϵ_max, 0, 1-εₚ)
            u = atanh(√a)
        end,
        u -> begin # de_du_derivative
            ∂ϵ╱∂u = 2ϵ_max * tanh(u) * sech(u)^2
        end
    )

end


"""Parametrization with a Logistic function that enforces ϵ_min < ϵ(t) < ϵ_max.
"""
function LogisticParametrization(ϵ_min, ϵ_max; k=1.0)

    Δ = ϵ_max - ϵ_min
    ε₀ = eps(0.0)  # 5⋅10⁻³²⁴
    @assert ϵ_max > ϵ_min

    PulseParametrization(
        u -> begin # epsilon_of_u
            ϵ = Δ / (1 + exp(-k * u)) + ϵ_min
        end,
        ϵ -> begin # u_of_epsilon
            ϵ′ = ϵ - ϵ_min
            a = max(ϵ′/(Δ - ϵ′), ε₀)
            u = log(a) / k
        end,
        u -> begin # de_du_derivative
            e⁻ᵏᵘ = exp(-k * u)
            ∂ϵ╱∂u = Δ * k * e⁻ᵏᵘ / (1 + e⁻ᵏᵘ)^2
        end
    )

end


"""Parametrization with a Logistic-Square function that enforces 0 ≤ ϵ(t) < ϵ_max.
"""
function LogisticSqParametrization(ϵ_max; k=1.0)

    ε₀ = eps(0.0)  # 5⋅10⁻³²⁴
    @assert ϵ_max > 0

    PulseParametrization(
        u -> begin # epsilon_of_u
            ϵ = ϵ_max * (2 / (1 + exp(-k * u)) - 1)^2
        end,
        ϵ -> begin # u_of_epsilon
            ρ = clamp(ϵ / ϵ_max, 0.0, 1.0)
            a = clamp((2 / (√ρ + 1)) - 1, ε₀, 1.0)
            u = -log(a) / k
        end,
        u -> begin # de_du_derivative
            eᵏᵘ = exp(k * u)
            ∂ϵ╱∂u = 4k * ϵ_max * eᵏᵘ * (eᵏᵘ  - 1) / (eᵏᵘ + 1)^3
        end
    )

end
