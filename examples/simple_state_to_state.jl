# # Optimization of a State-to-State Transfer in a Two-Level-System

#md # !!! tip
#md #     This example is also available as a Jupyter notebook:
#md #     [`simple_state_to_state.ipynb`](@__NBVIEWER_ROOT_URL__/examples/simple_state_to_state.ipynb)

# This first example illustrates the basic use of the `Krotov.jl` by solving a
# simple canonical optimization problem: the transfer of population in a two
# level system.

const σ̂_z = ComplexF64[1 0; 0 -1]
const σ̂_x = ComplexF64[0 1; 1  0]



"""Two-level-system Hamiltonian."""
function hamiltonian(Ω=1.0, E0=0.2)

    Ĥ₀ = -0.5 * Ω * σ̂_z
    Ĥ₁ = σ̂_x
    
    ϵ(t) = 0

    return [Ĥ₀, [Ĥ₁, ϵ]]

end
