# SEIR(S) epidemics on activity-structured population with time-varying individual behaviour

## Alexis Vargas Richards, Imperial College London Department of Infectious Disease Epidemiology

[![CI](https://github.com/vargasrichards/MixSwitchStochEpi.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/vargasrichards/MixSwitchStochEpi.jl/actions)
[![Docs build](https://github.com/vargasrichards/MixSwitchStochEpi.jl/actions/workflows/docs.yml/badge.svg)](https://github.com/vargasrichards/MixSwitchStochEpi.jl/actions)
[![Docs (GitHub Pages)](https://img.shields.io/badge/docs-GitHub%20Pages-blue.svg)](https://vargasrichards.github.io/MixSwitchStochEpi.jl)
[![Codecov](https://codecov.io/gh/vargasrichards/MixSwitchStochEpi.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/vargasrichards/MixSwitchStochEpi.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![License](https://img.shields.io/github/license/vargasrichards/MixSwitchStochEpi.jl)](LICENSE)
[![Status](https://img.shields.io/badge/status-WIP-orange.svg)](https://github.com/vargasrichards/MixSwitchStochEpi.jl)
[![Julia](https://img.shields.io/badge/julia-1.12-orange.svg)](https://julialang.org/)

## Overview

This repository implements a class-structured SEIRS epidemic model with both deterministic ODE and stochastic Gillespie-style SSA flavours. It includes utilities for constructing mixing and switching matrices, parameter calibration, and tools to run parameter sweeps and extinction/fadeout analyses.

The codebase follows a script-first workflow (files use `include(...)`), so examples below assume you run them from the repository root with the project activated.

## Installation

Ensure you have Julia (recommended 1.12). From the repository root run:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

## Quickstart — stochastic single run

Open a Julia REPL in the repo root and run:

```julia
using MixSwitchStochEpi

# Get default parameters + initial state
params, init, tspan = default_parametrisation()

# Run a single stochastic trajectory
res = simulate_system!(params, init; tmax = 100.0)

println("n snapshots: ", length(res.times))
display(res.states[end])  # final S,E,I,R per class
```

## Quickstart — many replicates & fadeout probability

```julia
using Statistics

results = many_simulations(params, init, 1, 200; tmax = 200.0)
println("mean secondary size: ", mean(results.secondary_sizes))

pext = pr_fadeout(params, init; n_sims = 1000, t_cutoff = 50.0)
println("Pr(extinction) ≈ ", pext)
```

Notes
- `simulate_system!` is the production stochastic simulator used throughout the tests.
- Deterministic ODE helpers reside in the repository and can be run via `include("src/classed_model.jl")` (see files in `src/`).

## Testing

Run the test suite:

```bash
julia --project=. test/runtests.jl
```

The test suite covers mixing matrices, switching, extinction behaviour, GPU vs CPU consistency, and other model invariants.

### CI and performance notes

- The full test-suite contains computationally heavy Monte Carlo and GPU tests. CI is configured to run a faster profile on PRs and the full suite on merges/nightly.
- To speed up PR runs locally or in CI use environment variables (CI_MODE, N_SIMS) — see test files for usage.

## Contributing

Please open issues or PRs. When adding/removing dependencies, update `Project.toml`/`Manifest.toml` and ensure `Pkg.test()` passes.

## License

This code is released under the MIT License — see `LICENSE`.

## Planned development

In the future, a spatially-embedded individual-based simulator will be added to investigate temporal activity variation and phylodynamic consequences.

## Citation

If you use this code in published work please cite the repository. A small BibTeX file is included: `julia_citations.bib`.