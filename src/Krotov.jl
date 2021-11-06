module Krotov

include("parametrization.jl")
export SquareParametrization, TanhParametrization, TanhSqParametrization, LogisticParametrization, LogisticSqParametrization

include("optimize.jl")
export optimize_pulses


import QuantumControlBase: optimize

"""
```julia
opt_result = optimize(problem; method=:krotov, kwargs...)
```

optimizes `problem` using Krotov's method, see
[`Krotov.optimize_pulses`](@ref).
"""
optimize(problem, method::Val{:krotov}; kwargs...) = Krotov.optimize_pulses(problem, kwargs...)


end
