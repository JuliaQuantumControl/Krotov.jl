module Krotov

include("parametrization.jl")
export SquareParametrization, TanhParametrization, TanhSqParametrization, LogisticParametrization, LogisticSqParametrization

include("result.jl")
include("workspace.jl")
include("optimize.jl")

import QuantumControlBase

"""
```julia
opt_result = optimize(problem; method=:krotov, kwargs...)
```

optimizes [`problem`](@ref QuantumControlBase.ControlProblem) using Krotov's
method, see [`Krotov.optimize_krotov`](@ref).
"""
QuantumControlBase.optimize(problem, method::Val{:krotov}; kwargs...) = optimize_krotov(problem, kwargs...)
QuantumControlBase.optimize(problem, method::Val{:Krotov}; kwargs...) = optimize_krotov(problem, kwargs...)


end
