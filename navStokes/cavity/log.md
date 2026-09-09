# navStokes/cavity — Log

_Started: 2026-09-04 | Last updated: 2026-09-09_

## Action Items

- [ ] Re=3200/5000 mismatch vs. Ghia is more likely **mesh
      under-resolution than a turbulence-modeling problem**. Ghia et
      al. (1982)'s benchmark is laminar at *every* Re they tabulate,
      including 5000 and 10000 — a strictly 2D lid-driven cavity stays
      steady/laminar well past Re=5000 (no vortex-stretching mechanism
      without a third dimension); the literature's Hopf bifurcation to
      unsteadiness is generally cited around Re≈8000–10000, and even
      that is periodic, not turbulent. Adding a RAS/turbulence model
      here would make the comparison invalid (RANS-averaged result vs.
      a laminar reference solution), the same reasoning that ruled out
      RANS in the old `channel_flow` project.
- [ ] `mesh/mesh.py` still uses a **uniform 50-point transfinite mesh,
      no wall clustering** (`points = 50`, no `Bump`/`Progression`) for
      every Re from 100 to 5000. As Re increases, near-wall shear
      layers thin and corner eddies shrink/multiply — Ghia's own
      higher-Re solutions needed a finer graded grid (129×129) plus a
      grid-independence check to resolve those features. Likely root
      cause of the growing Re=3200/5000 mismatch; consider a graded
      (wall-clustered) mesh if pursuing those Re further.
- [ ] Before refining the mesh, rule out non-convergence first (cheaper
      check): `residualControl` was only added to `fvSolution` on
      2026-09-09 — the existing Re=3200/5000 runs may predate it and
      simply not have run long enough to reach steady state (viscous
      diffusion time scales with Re, so higher-Re cases need
      proportionally more physical time to settle).
- [ ] `Allrun` still runs `blockMesh`+`mirrorMesh` (stale, copied from
      the now-deleted `cylinderLowRe` template) instead of the real
      `mesh.py`+`gmshToFoam` workflow.
- [ ] Consider retyping the `lid` physical surface to `wall` if
      wall-function/`yPlus` post-processing is ever needed.
- [ ] `plots.gp` output is raster PNG; vector (`pdfcairo`/`epslatex`)
      considered but not applied — see `plot_notes.md`.
- [ ] This session's `.gitignore` fix + reorg are still uncommitted —
      review before committing.

## Objective

OpenFOAM case for classic 2D lid-driven cavity flow — unit square, top
("lid") wall moving tangentially, other three walls stationary. Mesh
built with `mesh.py` (gmsh Python API, thin `d/100` single-layer
extrusion for OpenFOAM's pseudo-2D convention), converted via
`gmshToFoam`, run with `foamRun -solver incompressibleFluid`. Validated
against Ghia, Ghia & Shin (1982)'s tabulated centerline profiles across
Re=100/400/1000/3200/5000.

As of 2026-09-07 this is the only case under `navStokes/` — sibling cases
`channel_flow/` and `cylinderLowRe/` were deleted by the user (confirmed
intentional). Layout: `navStokes/cavity/` holds `mesh.py`, `log.md`,
`plots.gp` (+ `plot_notes.md`), `benchmarks/` (Ghia data), `plots/`
(figures); the actual OpenFOAM case lives in `navStokes/cavity/openFoam/`
(was `navStokes/openFoam/cavity/` before the 2026-09-07 reorg).

## Key decisions

- **Mesh/case setup**: `Allclean` is bash (container runs bash, not
  fish). Physical groups: `lid`/`wall`/`empty`(front-back caps)/`domain`,
  matching `channel_flow/mesh.geo`'s old naming convention.
- **`oforg.fish`** now runs the container as `--user (id -u):(id -g) -e
  HOME=/tmp` — bind-mounted output was coming out root-owned otherwise.
- **`.gitignore`** generalized `navStokes/openFoam/**/...` →
  `navStokes/**/openFoam/**/...` so the ignore rules survive `openFoam`
  moving from container-directory to case-subdirectory (broke silently
  after the 2026-09-07 reorg; ~653 generated files had gone from ignored
  to staged before the fix).
- **2026-09-09** — Added `residualControl { U 1e-5; p 1e-5; }` to
  `system/fvSolution`'s `PIMPLE` block, to auto-stop the run once the
  solution settles to steady state instead of always running to
  `endTime`.
  Why: `solvers{}`'s `tolerance`/`relTol` only control how tightly each
  linear system is solved *within* one time step (numerical accuracy per
  step) — they don't detect whether the physical solution has stopped
  changing between steps. `residualControl` is the actual mechanism for
  that: it checks each listed field's *initial* residual every time step
  and calls an early stop once all of them are below threshold
  simultaneously.

## Learnings (technical gotchas worth not re-deriving)

- `gmsh.model.geo.extrude(...)` (Python API) returns `[far cap, volume,
  lateral surfaces in Curve Loop order]` — same convention as `.geo`
  `Extrude{}`. Verified against `constant/polyMesh/boundary` face counts.
- `gmshToFoam` infers boundary `type` (`wall` vs. generic `patch`) from
  the physical-surface *name* — `"lid"` got `type patch` despite being a
  real wall. Harmless for laminar solving; would matter for
  wall-function/`yPlus` post-processing.
- OpenFOAM's `internalField` is per-*cell* (mesh-numbering order), not a
  reshape-able i,j grid. Use the `sample` function object
  (`system/sampleDict`, `type sets;`, `postProcess -func sampleDict`) for
  `x=const`/`y=const` profiles — it interpolates via the mesh's own
  connectivity.
- Validation-plot convention: tabulated benchmark data as discrete
  markers, current simulation as a line — not both `with line`. With
  multiple Re on one figure, differentiate further within each family:
  marker shape per Re for Ghia, dash pattern per Re for OpenFOAM.
- `git status`'s `R` (rename) entries can pair an unrelated staged
  deletion + addition purely by content similarity — not proof of an
  actual move. Untracked files also don't travel with a `git`-aware
  directory reorg (left `plot_notes.md`/`plots.gp.bak` stranded at the
  old path; `plots.gp.bak` ended up in the OS Trash — recovered).
- `.gitignore` prefixes anchored to a literal path break silently on
  rename; anchor on a path *segment name* (`**/openFoam/**`) instead.
- `tolerance`/`relTol` (in `fvSolution`'s `solvers{}`) vs.
  `residualControl` (in `PIMPLE{}`/`SIMPLE{}`) are not interchangeable:
  the former is per-time-step linear-solver accuracy, the latter is
  physical-convergence/steady-state detection across time steps. Easy to
  conflate since both look like "convergence tolerance" settings.
- Ghia et al. (1982) is laminar at every tabulated Re (100–10000) — 2D
  lid-cavity flow doesn't turbulence-transition the way 3D flow does;
  don't reach for a turbulence model to fix a high-Re mismatch against
  this benchmark specifically.

## References

- Ghia, U., Ghia, K.N., Shin, C.T. (1982), *High-Re solutions for
  incompressible flow using the Navier-Stokes equations and a multigrid
  method*, J. Comput. Phys. 48(3):387–411.
