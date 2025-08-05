# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Krotov.jl is a Julia package implementing Krotov's method of optimal control for quantum systems, enhanced with automatic differentiation. It is part of the JuliaQuantumControl organization and designed to work with the QuantumControl.jl framework.

The package is a port of the krotov Python package, adapted to the API of QuantumControl.jl.

## Development Commands

### Development Environment
- `make devrepl`: Start an interactive REPL for testing and building documentation (recommended)
- `julia -i --project=test devrepl.jl`: Alternative way to start development REPL

### Testing
- `make test`: Run the complete test suite
- `julia --project=test --banner=no --startup-file=yes -e 'include("devrepl.jl"); test()'`: Alternative test command
- Test individual files by running them from the test REPL

### Documentation
- `make docs`: Build the documentation

### Code Quality
- `make codestyle`: Apply JuliaFormatter to the entire project
- Requires `../.JuliaFormatter.toml` configuration file

### Cleanup
- `make clean`: Clean up build/doc/testing artifacts
- `make distclean`: Restore to a clean checkout state

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

- Part of JuliaQuantumControl ecosystem - may use shared development scripts in `../scripts/`
- This package is designed to work within the JuliaQuantumControl development environment
- Code formatting follows JuliaQuantumControl organization standards
- Tests require the full test environment with additional dependencies
- Uses `devrepl.jl` for unified development environment setup
