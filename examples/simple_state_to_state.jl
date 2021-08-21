# # Optimization of a State-to-State Transfer in a Two-Level-System

#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`simple_state_to_state.ipynb`](@__NBVIEWER_ROOT_URL__/examples/simple_state_to_state.ipynb)

# This first example illustrates the basic use of the `Krotov.jl` by solving a
# simple canonical optimization problem: the transfer of population in a two
# level system.

using QuantumControlBase
using Krotov
using Plots
#-
const σ̂_z = ComplexF64[1 0; 0 -1]
const σ̂_x = ComplexF64[0 1; 1  0]


#-
"""Two-level-system Hamiltonian."""
function hamiltonian(Ω=1.0, E0=0.2)

    Ĥ₀ = -0.5 * Ω * σ̂_z
    Ĥ₁ = σ̂_x

    ϵ(t) = E0 * flattop(t, T=5, t_rise=0.3, func=:blackman)

    return [Ĥ₀, [Ĥ₁, ϵ]]

end
#-

H = hamiltonian()
#jl @test length(H) == 2
#-
tlist = collect(range(0, 5, length=500))
#-

plot(tlist, H[2][2].(tlist))


# The control field here switches on from zero at $t=0$ to it's maximum amplitude
# 0.2 within the time period 0.3 (the switch-on shape is half a [Blackman pulse](https://en.wikipedia.org/wiki/Window_function#Blackman_window)).
# It switches off again in the time period 0.3 before the
# final time $T=5$). We use a time grid with 500 time steps between 0 and $T$:
