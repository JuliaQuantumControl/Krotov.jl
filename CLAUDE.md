# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Krotov.jl is a Julia package implementing Krotov's method of optimal control for quantum systems, enhanced with automatic differentiation. It is part of the JuliaQuantumControl organization and designed to work with the QuantumControl.jl framework.

The package is a port of the krotov Python package, adapted to the API of QuantumControl.jl.

## Development Commands

Run `make help` for all targets. The development workflow is documented in the org-wide [CONTRIBUTING.md](https://github.com/JuliaQuantumControl/.github/blob/master/CONTRIBUTING.md) (`../.github/CONTRIBUTING.md` in the development environment).

- `make test`: Run the test suite in the `test` environment (or `julia --project=test -e 'include("test/runtests.jl")'`)
- `make devrepl`: REPL with the `test` environment active and the `docs` environment stacked; run individual test files (`include("test/test_tls_optimization.jl")`), `include("test/runtests.jl")`, or `include("docs/make.jl")` from there
- `make docs`: Build the documentation in the `docs` environment
- `make coverage` / `make htmlcoverage`: Test coverage
- `make codestyle`: Apply JuliaFormatter (version pinned in the `Makefile`) and check `[sources]`
- `make clean` / `make distclean`: Remove build/test artifacts

Sibling packages (QuantumControl, QuantumPropagators, GRAPE, …) come from their registered releases, or temporarily from a GitHub branch via a URL `[sources]` entry in `test/Project.toml` / `docs/Project.toml`. Never commit a `path` source for a sibling (as written by `../scripts/installorg.jl`). The `test` and `docs` environments reference the package itself via `[sources]` (`{path = ".."}`); this needs Julia ≥ 1.11.

## Architecture

### Core Module Structure
- `src/Krotov.jl`: Main module file that includes other components
- `src/optimize.jl`: Main optimization logic implementing Krotov's method
- `src/result.jl`: `KrotovResult` type for storing optimization results
- `src/workspace.jl`: `KrotovWrk` workspace for internal optimization state

### Key Components

#### Optimization Flow
The optimization follows Krotov's method with forward/backward propagation:
1. Forward propagation of initial states
2. Backward propagation with adjoint states
3. Pulse updates using gradient information
4. Iteration until convergence

#### Result Object (`KrotovResult`)
Stores optimization results including:
- Iteration information and convergence status
- Final-time functional values (`J_T`, `J_T_prev`)
- Original and optimized control fields
- Target state overlaps (`tau_vals`)
- Timing and callback records

#### Workspace (`KrotovWrk`)
Internal workspace containing:
- Trajectories and adjoint trajectories
- Forward/backward propagators and storage
- Control derivatives and pulse discretizations
- Lambda values and update shapes

### Integration Points
- Extends `QuantumControl.optimize` with `method=Krotov`
- Uses `QuantumControl.QuantumPropagators` for state propagation
- Implements `AbstractOptimizationResult` interface
- Supports threading via `@threadsif` for parallel trajectory propagation

### Test Structure
- `test/runtests.jl`: Main test runner using SafeTestsets
- Individual test files for different optimization scenarios:
  - `test_tls_optimization.jl`: Two-level system optimization
  - `test_pulse_optimization.jl`: General pulse optimization
  - `test_iterations.jl`: Iteration mechanics
  - `test_empty_optimization.jl`: Edge cases

## Development Notes

- Part of the JuliaQuantumControl ecosystem
- Code formatting follows JuliaQuantumControl organization standards
- Tests use `QuantumControl.DummyOptimization` (experimental) for dummy control problems and `QuantumControlTestUtils.RandomObjects` for random states and matrices
