```@meta
CurrentModule = Krotov
```

# Krotov.jl

```@eval
using Markdown
using Pkg

VERSION = Pkg.dependencies()[Base.UUID("b05dcdc7-62f6-4360-bf2c-0898bba419de")].version

github_badge = "[![Github](https://img.shields.io/badge/JuliaQuantumControl-Krotov.jl-blue.svg?logo=github)](https://github.com/JuliaQuantumControl/Krotov.jl)"

version_badge = "![v$VERSION](https://img.shields.io/badge/version-v$VERSION-green.svg)"

Markdown.parse("$github_badge $version_badge")
```

Implementation of [Krotov's method of optimal control](https://arxiv.org/abs/1008.5126) enhanced with automatic differentiation.

Part of [`QuantumControl.jl`](https://github.com/JuliaQuantumControl/QuantumControl.jl#readme) and the [JuliaQuantumControl](https://github.com/JuliaQuantumControl) organization.

[Krotov.jl](https://github.com/JuliaQuantumControl/Krotov.jl) is a port of the [`krotov` Python package](https://github.com/qucontrol/krotov#readme), adapted to the API  of [`QuantumControl.jl`](https://github.com/JuliaQuantumControl/QuantumControl.jl#readme).


## Contents

### Overview

```@contents
Pages = [
    "overview.md",
]
Depth = 1
```

### Examples

```@contents
Pages = [
    "examples/simple_state_to_state.md",
    "examples/rho_3states.md",
    "examples/state_to_state_parametrizations.md",
    "examples/perfect_entanglers.md",
]
Depth = 1
```

See also the [general examples](https://juliaquantumcontrol.github.io/QuantumControl.jl/stable/examples/) of the [QuantumControl](https://juliaquantumcontrol.github.io/QuantumControl.jl/stable/) package.


### API

```@contents
Pages = [
    "api.md",
]
Depth = 1
```

## History

See the [Releases](https://github.com/JuliaQuantumControl/Krotov.jl/releases) on Github.

## References

```@bibliography
```
