# Quick orientation for AI coding agents

This repository implements SEIR-style epidemic models in R
(package-style) with supplemental Julia code for stochastic simulations.
The guidance below is concise and immediately actionable for making code
edits, tests, docs, and small feature additions.

## Big picture

- R package live under `R/` with documentation in `man/` and pkgdown
  site built into `docs/`.
- Tests use `testthat/` (see
  `tests/testthat/test-seirs_mixing_switching.R`).
- Julia stochastic models and visualization code live under `julia/`
  (see `julia/stochastic_classed/SEIRS_stochastic.jl`).
- Dependencies: R uses `renv` (see `renv.lock`) and the package `odin2`
  is required; Julia deps are in `Project.toml` under the repo root
  (used for `julia/` scripts).

## Key files and examples to reference

- `R/parameter_store.R` — central place for scenario/parameter handling.
- `R/seir_mixswitch.R` (and `R/unstructured_seirs.R`) — deterministic
  model implementations.
- `tests/testthat/test-seirs_mixing_switching.R` — canonical test
  showing expected outputs and shapes.
- `julia/stochastic_classed/SEIRS_stochastic.jl` — stochastic simulation
  implementation; follow its data shapes when producing simulated
  outputs.
- `docs/llms.txt` and `README.md` — short human-readable repo
  descriptions to reuse or quote.

## Developer workflows (explicit, reproducible steps)

- Recreate R environment: open R and run
  [`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html)
  to install packages pinned in `renv.lock`.
- Run R package tests: in R,
  [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
  or `testthat::test_dir('tests/testthat')`.
- Regenerate documentation after code changes:
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
  (updates NAMESPACE and `man/`).
- Build the website/docs:
  [`pkgdown::build_site()`](https://pkgdown.r-lib.org/reference/build_site.html)
  will update `docs/`.
- Run R CMD check for package-level validation: from shell
  `R CMD check .` (or
  [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
  in R).
- Julia environment: run
  `julia --project=julia -e 'using Pkg; Pkg.instantiate()'` to install
  Julia deps listed, then execute Julia scripts with
  `julia --project=julia path/to/script.jl`.

## Project-specific conventions

- Documentation: use roxygen2-style comments in `R/` and call
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html);
  doc files in `man/` are used to feed pkgdown.
- Tests: follow `testthat` structure; tests validate numerical outputs
  and object shapes rather than extensive edge-case coverage.
- Naming: many functions refer to “activity classes” and “switching” —
  preserve the term usage (e.g. `mean_node_switching`,
  `visualise_mixswitch`).
- Output files and plots are kept at repo root (e.g. `*.pdf`, `*.svg`)
  and in `docs/` — prefer not to commit large binary outputs when
  editing code unless intentionally updating figures.

## Integration points and external dependencies

- R: `odin2` (runtime for some model components) — ensure it’s installed
  in the R environment.
- Julia: packages listed in `Project.toml` (Agents, Catalyst, Makie
  stacks) — used by `julia/` scripts.
- Data/parameters: look at `parameters.csv`, `output.csv` at repo root
  for example formats.

## Practical guidance for AI code edits

- When adding or changing exported R functions:
  - Edit file in `R/`, add roxygen comments, run
    [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
    to update `man/` and `NAMESPACE`.
  - Add or update a focused `testthat` test in `tests/testthat/`
    mirroring existing test style (`test-seirs_mixing_switching.R`).
- When touching Julia scripts, prefer small, self-contained changes and
  keep the same project environment (`--project=julia`).
- For changes that affect docs or examples, run
  [`pkgdown::build_site()`](https://pkgdown.r-lib.org/reference/build_site.html)
  and include only the updated HTML/MD files in `docs/` if the change is
  intended for publication.

## What to avoid / assumptions

- Do not assume a global R or Julia environment; always document
  environment changes and prefer using
  [`renv::restore()`](https://rstudio.github.io/renv/reference/restore.html)
  and `Pkg.instantiate()`.
- Avoid adding heavy binary blobs; use the existing pattern of generated
  plots committed only when they are intentional deliverables.

If any part of the runtime setup is unclear (specific odin2 models used,
or a Julia script’s expected CLI inputs), tell me which file you want to
modify and I will extract exact runtime calls and add a short runnable
example to this instruction file.
