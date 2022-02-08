# Krotov.jl

[![Version](https://juliahub.com/docs/Krotov/version.svg)](https://juliahub.com/ui/Packages/Krotov/3mCxK)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaquantumcontrol.github.io/Krotov.jl/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaquantumcontrol.github.io/Krotov.jl/dev)
[![Build Status](https://github.com/JuliaQuantumControl/Krotov.jl/workflows/CI/badge.svg)](https://github.com/JuliaQuantumControl/Krotov.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaQuantumControl/Krotov.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaQuantumControl/Krotov.jl)

Implementation of [Krotov's method of optimal control](https://arxiv.org/abs/1008.5126), part of [`QuantumControl.jl`][QuantumControl] and the [JuliaQuantumControl][] organization.

This package is a port of the [`krotov` Python package](https://github.com/qucontrol/krotov#readme), adapted to the API  of [`QuantumControl.jl`][QuantumControl].

## Installation

For normal usage, the `Krotov` package should be installed as part of [`QuantumControl.jl`][QuantumControl]:

~~~
pkg> add QuantumControl
~~~

For development usage, see the [organization development notes](https://github.com/JuliaQuantumControl#development).

## Documentation

A minimal standalone documentation of `Krotov.jl` is available at <https://juliaquantumcontrol.github.io/Krotov.jl>.

For a broader perspective, see the [documentation of the `QuantumControl.jl` meta-package](https://juliaquantumcontrol.github.io/QuantumControl.jl/).

See also the [documentation of the `krotov` Python package](https://qucontrol.github.io/krotov) for an in-depth overview of the core method.

[QuantumControl]: https://github.com/JuliaQuantumControl/QuantumControl.jl#readme
[JuliaQuantumControl]: https://github.com/JuliaQuantumControl
