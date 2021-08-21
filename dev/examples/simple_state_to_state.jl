using QuantumControlBase
using Krotov
using PyPlot

const σ̂_z = ComplexF64[1 0; 0 -1]
const σ̂_x = ComplexF64[0 1; 1  0]

"""Two-level-system Hamiltonian."""
function hamiltonian(Ω=1.0, E0=0.2)

    Ĥ₀ = -0.5 * Ω * σ̂_z
    Ĥ₁ = σ̂_x

    ϵ(t) = E0 * flattop(t, T=5, t_rise=0.3, func=:blackman)

    return [Ĥ₀, [Ĥ₁, ϵ]]

end

H = hamiltonian();
@test length(H) == 2

tlist = collect(range(0, 5, length=500));

plot(tlist, H[2][2].(tlist));

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

