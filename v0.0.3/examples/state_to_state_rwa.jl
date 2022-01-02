const σ̂_z = ComplexF64[1 0; 0 -1]
const σ̂_x = ComplexF64[0 1; 1  0]


"""Two-level-system Hamiltonian."""
function hamiltonian(Ω=1.0, E0=0.2)

    Ĥ₀ = -0.5 * Ω * σ̂_z
    Ĥ₁ = σ̂_x

    return [Ĥ₀, Ĥ₁]

end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

