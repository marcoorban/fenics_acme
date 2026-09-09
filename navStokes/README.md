# navStokes

FEM solver for the incompressible Navier-Stokes equations, built in
FEniCS (`dolfinx`), with boundary conditions enforced **weakly** (Nitsche
method) rather than as strong Dirichlet constraints.

## Current scope

For now, the project targets two classic, well-documented benchmark
flows, each run under a range of flow conditions:

- **Lid-driven cavity flow** — see `cavity/`.
- **Flow past a cylinder** — at various Reynolds numbers.

Both are chosen because they have established baselines to validate
against: an OpenFOAM run of the same case (see `cavity/openFoam/`) plus
published benchmark data (e.g. Ghia, Ghia & Shin (1982) for the cavity —
`cavity/benchmarks/`). The FEM results are checked against both.

## Deferred: turbulent cavity flow

Turbulent (DNS-scale) lid-driven cavity flow is a planned extension, not
current work — it needs mesh resolution and compute far beyond the
laminar benchmark cases above, and isn't practical to take on yet.
