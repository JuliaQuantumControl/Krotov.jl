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

Implementation of Krotov's method of optimal control [Krotov1996,SomloiCP1993,BartanaJCP1997,PalaoPRA2003,ReichJCP2012,GoerzSPP2019](@cite) enhanced with automatic differentiation [GoerzQ2022](@cite).

Part of [`QuantumControl.jl`](https://github.com/JuliaQuantumControl/QuantumControl.jl#readme) and the [JuliaQuantumControl](https://github.com/JuliaQuantumControl) organization.

[Krotov.jl](https://github.com/JuliaQuantumControl/Krotov.jl) is a port of the [`krotov` Python package](https://github.com/qucontrol/krotov#readme), adapted to the API  of [`QuantumControl.jl`](https://github.com/JuliaQuantumControl/QuantumControl.jl#readme).


## Contents

```@contents
Depth = 2
Pages = [pair[2] for pair in Main.PAGES[2:end-1]]
```


## History

See the [Releases](https://github.com/JuliaQuantumControl/Krotov.jl/releases) on Github.
