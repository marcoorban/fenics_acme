# navStokes/cavity — Log

_Started: 2026-09-04 | Last updated: 2026-09-10_

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
- [ ] 2026-09-10: `mesh/mesh.py` bumped from 50 to **130 points
      (129×129 cells)** — now matches Ghia's own grid resolution at
      their higher Re. Still **uniform, no wall clustering**
      (`Bump`/`Progression`); as Re increases, near-wall shear layers
      thin and corner eddies shrink further, and Ghia's own solutions
      still relied on grid-independence checks even at 129×129, so this
      may not fully resolve the Re=3200/5000 mismatch on its own — worth
      re-plotting against Ghia at this resolution before deciding
      whether wall clustering is still needed.
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
- [ ] 2026-09-10 parallel-run test left real (non-throwaway) output:
      `processor0-3/` and a converged `t≈1.7334` result from whatever
      Re was configured at the time, produced while verifying the
      `decomposePar`/`mpirun`/`reconstructPar` pipeline end-to-end. Not
      cleaned up — decide whether to keep or clear it.

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
- **2026-09-10** — Case can now run in parallel: added
  `system/decomposeParDict` and updated `foamRun.fish` to
  `decomposePar && mpirun -np N foamRun -solver incompressibleFluid
  -parallel && reconstructPar -latestTime && foamPostProcess -func
  sampleDict -latestTime`.
  Why 4 initially, then 8: machine has 20 logical threads (i7-12700,
  12 cores × 2-way hyperthreading). Started at 4 (conservative,
  matched the then-2401-cell mesh); bumped to 8
  (`numberOfSubdomains 8; simpleCoeffs { n (4 2 1); ... }` — a 4×2
  split in x/y, 1 in z) the same day the mesh grew to 129×129
  (16641 cells), which justifies more subdomains. `-np` here must
  always match `numberOfSubdomains` if either changes.
  Verified end-to-end by actually running it both times, not just
  reading the dictionary syntax: at `-np 4`, `decomposePar` → parallel
  solve (converged via `residualControl` at t≈1.73s) → `reconstructPar`
  → `foamPostProcess` all completed with exit code 0; at `-np 8` on the
  129×129 mesh, `decomposePar` alone was re-verified (subdomains came
  out well balanced: 2079–2081 cells each, despite 129 not dividing
  evenly by 4 or 2 — `simple` decomposition absorbs the remainder
  gracefully).

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
  multiple Re on one figure, differentiate further with one color per
  Re (Okabe-Ito colorblind-safe palette), shared between a Ghia marker
  and its matching OpenFOAM line, plus marker shape per Re as a
  redundant cue. Dash patterns per Re were tried first and dropped —
  with 5 overlapping curves, distinguishing dash patterns at crossings
  was harder to read than just using color.
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
- `decomposePar` decomposes *every* field file present in `0/`,
  regardless of whether the active `simulationType` actually reads it.
  `0/nut`, `0/k`, `0/epsilon`, `0/omega`, `0/nuTilda` were leftover
  stock-tutorial files (wrong patch names too: `movingWall`/
  `fixedWalls`/`frontAndBack` instead of this case's `lid`/`wall`/
  `empty`) — harmless in a serial `simulationType laminar` run since
  the solver never touches them, but fatal to `decomposePar`. Deleted
  all five (2026-09-10; same cleanup already done for `channel_flow`).
  Worth checking `0/` for this kind of dead tutorial leftover on any
  future case before trying to parallelize it.

## References

- Ghia, U., Ghia, K.N., Shin, C.T. (1982), *High-Re solutions for
  incompressible flow using the Navier-Stokes equations and a multigrid
  method*, J. Comput. Phys. 48(3):387–411.
