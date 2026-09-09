# navStokes/cavity — Log

_Started: 2026-09-04 | Last updated: 2026-09-07 (condensed)_

## Objective

OpenFOAM case for classic 2D lid-driven cavity flow — unit square, top
("lid") wall moving tangentially, other three walls stationary. Mesh
built with `mesh.py` (gmsh Python API, thin `d/100` single-layer
extrusion for OpenFOAM's pseudo-2D convention), converted via
`gmshToFoam`, run with `foamRun -solver incompressibleFluid`. Validated
against Ghia, Ghia & Shin (1982)'s tabulated centerline profiles.

Currently Re=100 (`constant/physicalProperties` `nu=1e-2`), `deltaT=0.01`
(Co≈0.49) — laminar, matches Ghia's benchmark table directly.

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
  markers, current simulation as a line — not both `with line`.
- `git status`'s `R` (rename) entries can pair an unrelated staged
  deletion + addition purely by content similarity — not proof of an
  actual move. Untracked files also don't travel with a `git`-aware
  directory reorg (left `plot_notes.md`/`plots.gp.bak` stranded at the
  old path; `plots.gp.bak` ended up in the OS Trash — recovered).
- `.gitignore` prefixes anchored to a literal path break silently on
  rename; anchor on a path *segment name* (`**/openFoam/**`) instead.

## Future Action Items

- [ ] `Allrun` still runs `blockMesh`+`mirrorMesh` (stale, copied from
      the now-deleted `cylinderLowRe` template) instead of the real
      `mesh.py`+`gmshToFoam` workflow.
- [ ] Confirm whether the run stopping at `t=161` (of `endTime 200`) is
      intentional/converged before treating `postProcessing/sampleDict/
      161/` as final.
- [ ] Consider retyping the `lid` physical surface to `wall` if
      wall-function/`yPlus` post-processing is ever needed.
- [ ] `plots.gp` output is raster PNG; vector (`pdfcairo`/`epslatex`)
      considered but not applied — see `plot_notes.md`.
- [ ] This session's `.gitignore` fix + reorg are still uncommitted —
      review before committing.

## References

- Ghia, U., Ghia, K.N., Shin, C.T. (1982), *High-Re solutions for
  incompressible flow using the Navier-Stokes equations and a multigrid
  method*, J. Comput. Phys. 48(3):387–411.
