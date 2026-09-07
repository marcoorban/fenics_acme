# navStokes/cavity — Log

_Started: 2026-09-04 | Last updated: 2026-09-07_

## Status

- **2026-09-04** — Re and `deltaT` findings below (flagged, left
  unchanged at the time) have since been fixed by the user: `nu` is now
  `1e-02` (Re=100, standard Ghia et al. laminar benchmark case) and
  `deltaT` is `0.01` (Courant ≈0.49, comfortably under CFL≈1). Re-checked
  and confirmed — see Learnings for the numbers this superseded.

## Objective

OpenFOAM case for classic 2D lid-driven cavity flow — a unit square, top
("lid") wall moving tangentially at a fixed speed, other three walls
stationary. Mesh is built with a custom `mesh.py` (gmsh Python API) rather
than `blockMesh` or a `.geo` script, converted to OpenFOAM's
`constant/polyMesh` via `gmshToFoam`, and run with `foamRun -solver
incompressibleFluid`. Validated against Ghia, Ghia & Shin (1982)'s
tabulated centerline velocity profiles.

As of 2026-09-07, this is the only case left under `navStokes/` — the
`channel_flow/` and `cylinderLowRe/` cases (formerly siblings under
`navStokes/openFoam/`) were deleted in this session's directory reorg,
confirmed intentional by the user (not reorg fallout). Project layout is
now `navStokes/cavity/` (this directory) holding `mesh.py`, `log.md`,
`plots.gp`, `benchmarks/` (Ghia reference data), and `plots/` (rendered
figures) directly, with the actual OpenFOAM dictionaries/case living one
level down in `navStokes/cavity/openFoam/` — previously it was
`navStokes/openFoam/cavity/`.

## Decisions

- **2026-09-04** — `Allclean` is a **bash** script, not fish, even though
  fish is used interactively elsewhere in this repo.
  Why: the `microfluidica/openfoam:org` Docker image runs bash, not fish —
  a script meant to execute inside/alongside OpenFOAM tooling needs to
  match that shell.
  Ruled out: a fish version (first draft) — logic was equivalent (remove
  all time directories except `0`) but wouldn't run in the container.

- **2026-09-04** — `mesh.py` extrudes the 2D square one unit in `z`
  (later changed by the user to `d/100`, a thin single-cell-thick slab)
  rather than staying a bare 2D surface.
  Why: OpenFOAM requires a 3D mesh even for a physically-2D problem; a
  thin single-layer extrusion is the standard convention for this.

- **2026-09-04** — Front/back cap surfaces are tagged as a separate
  physical group named `"empty"`, not lumped into `"wall"`.
  Why: this is a 2D flow problem — OpenFOAM's convention for pseudo-2D
  cases is an `empty` patch type on the direction that isn't solved,
  which requires `0/U` and `0/p` to set `type empty` on exactly those
  faces.

- **2026-09-04** — Physical groups on the mesh: `"lid"` (top face),
  `"wall"` (bottom/left/right), `"empty"` (front/back caps), `"domain"`
  (the volume) — matching the naming pattern already used in
  `../channel_flow/mesh.geo` (`Physical Surface`/`Physical Volume` with
  descriptive names, plus a volume group so `gmshToFoam` has an
  unambiguous cell region).

- **2026-09-04** — `~/.config/fish/functions/oforg.fish` (the Docker
  wrapper used to run OpenFOAM) now passes `--user (id -u):(id -g) -e
  HOME=/tmp`.
  Why: without `--user`, the container runs as root; since `-v
  $PWD:/project` is a bind mount sharing the host's user namespace,
  everything OpenFOAM wrote (time directories, mesh output) came out
  root-owned on the host, blocking normal edits/cleanup without `sudo`.
  `HOME=/tmp` is needed alongside it because the host UID has no entry in
  the container's `/etc/passwd`, so it has no real `$HOME` by default.
  Ruled out: chown-after-the-fact wrapper, or Docker `userns-remap` at
  the daemon level — both work but are more invasive than just running
  the container as the host UID directly.

- **2026-09-07** — Project restructured from `navStokes/openFoam/cavity/`
  (case dictionaries direct, no separation from analysis files) to
  `navStokes/cavity/` with the OpenFOAM case itself moved down into a
  `openFoam/` subdirectory, alongside `mesh.py`/`log.md`/`plots.gp`/
  `benchmarks/`/`plots/` as siblings rather than case-internal files.
  `navStokes/openFoam/README.md`, `pyproject.toml`, and `uv.lock` moved
  up to `navStokes/`. Done by the user directly (not scripted this
  session); effect is that analysis/plotting artifacts are no longer
  nested inside the OpenFOAM case tree itself.
  Ruled out (by omission, confirmed intentional): keeping `channel_flow/`
  and `cylinderLowRe/` — both were deleted outright as part of this
  reorg, not moved anywhere.

- **2026-09-07** — `.gitignore`'s OpenFOAM section generalized from
  `navStokes/openFoam/**/...` to `navStokes/**/openFoam/**/...`.
  Why: the reorg above moved the case from `navStokes/openFoam/cavity/`
  to `navStokes/cavity/openFoam/` — the old pattern was anchored to
  `openFoam` being the top-level container right under `navStokes/`, so
  it silently stopped matching anything once `openFoam` became a
  sub-directory of the case instead. Caught because `git status` showed
  ~653 files staged as new: all 161 time-step directories,
  `constant/polyMesh/*`, `postProcessing/`, and `of.foam` had gone from
  ignored to tracked. Generalizing to match an `openFoam` path segment
  at *any* depth under `navStokes/` survives this kind of reorg without
  needing to remember to update `.gitignore` by hand each time.
  Ruled out: leaving the old hardcoded path and just re-adding a second
  copy of the same rules for the new path — works today but breaks again
  identically on the next rename; the `**`-generalized version doesn't.

## Tasks

- [x] `Allclean` — removes all time directories except `0` (2026-09-04) —
      written in fish first, then rewritten in bash per the Decision
      above.
- [x] `mesh.py` — one-unit (later `d/100`) extrusion + physical group
      tagging (`lid`/`wall`/`empty`/`domain`) (2026-09-04)
- [x] `oforg.fish` — run as host user instead of root, verified by the
      user that newly-written files are no longer root-owned (2026-09-04)
- [x] Reviewed `0/U`, `0/p`, `system/controlDict` (plus `fvSolution`,
      `fvSchemes`, `constant/physicalProperties`,
      `constant/polyMesh/boundary`) for lid-driven-cavity correctness
      (2026-09-04) — see Learnings below for findings.
- [x] Re and `deltaT` re-checked after the user's own fix (2026-09-04,
      same day, before the fixes below): `nu` → `1e-02` (Re=100,
      standard Ghia benchmark case), `deltaT` → `0.01` (Co≈0.49). Both
      confirmed resolved — see Status above.
- [x] Case run to `t=161` (of `endTime 200`); `postProcessing/sampleDict/
      161/` populated. Added `system/sampleDict` (`type sets;`, two
      `uniform` line samples — vertical centerline at `x=0.5` for `u(y)`,
      horizontal centerline at `y=0.5` for `v(x)`, `z=0.005` mid-plane) to
      extract profiles comparable to Ghia et al. (1982) Tables I/II, run
      via `postProcess -func sampleDict -time 161` (2026-09-07).
- [x] `plots.gp` backed up (`plots.gp.bak`) and restyled for a
      publication-appropriate look — Okabe-Ito colorblind-safe palette,
      Ghia data as discrete open-circle markers vs. OpenFOAM as a solid
      line (was: both `with line`), italicized axis labels, minor tics,
      light gridlines, `Re = 100` label on-figure, fixed a wrong "Ghia
      (1981)" citation year and an inconsistent legend label. Rendered
      and visually confirmed both figures — OpenFOAM tracks the Ghia
      benchmark closely. Full rationale and options *not* applied
      (vector/PDF output, no-gridline variant, etc.) written up in
      `plot_notes.md` for review (2026-09-07).
- [x] Read the reorganized directory tree after the user's move; found
      and fixed the `.gitignore` staleness (see Decisions) and unstaged
      the ~653 erroneously-tracked generated files, leaving only
      hand-authored case files (`Allrun`, `Allclean`, `0/*`,
      `constant/momentumTransport`, `constant/physicalProperties`,
      `system/*`) staged (2026-09-07).
- [x] Recovered `plot_notes.md` and `plots.gp.bak` — both were untracked
      when the reorg ran, so a `git`-based move left them stranded at the
      old `navStokes/openFoam/cavity/` path (and `plots.gp.bak`
      specifically ended up in the OS Trash). Moved both into
      `navStokes/cavity/`; old now-empty `navStokes/openFoam/` directory
      removed (2026-09-07).

## Learnings

- **2026-09-04** — `gmsh.model.geo.extrude([(2, surface)], dx, dy, dz,
  ...)` (Python API) returns entities in the same order as the `.geo`
  `Extrude{}` command already documented in `../channel_flow/log.md`:
  `out[0]` = far cap, `out[1]` = volume, `out[2..]` = lateral surfaces in
  `Curve Loop` order (bottom, right, top, left for this geometry).
  Confirmed by generating `cavity.msh` and checking both its
  `$PhysicalNames` block and the resulting `constant/polyMesh/boundary`
  face counts (`empty`=4802 = 2×49×49, `wall`=147 = 3×49, `lid`=49 —
  all match the expected face counts for a 50-point transfinite edge).

- **2026-09-04** — `gmshToFoam` infers the OpenFOAM boundary `type`
  (`wall` vs. generic `patch`) heuristically from the physical-surface
  *name* — a group literally named `"wall"` gets `type wall` in
  `constant/polyMesh/boundary`, but `"lid"` (despite being a real, moving
  wall) gets `type patch`. Numerically harmless for a laminar solver
  (the Dirichlet velocity BC behaves identically either way), but would
  matter if wall-function, `wallShearStress`, or `yPlus`
  post-processing is ever needed on the lid.

- **2026-09-04** — Docker bind mounts (`-v $PWD:/project`) share the
  host's user namespace: a container run without `--user` executes as
  root, so anything it writes to the mount comes out root-owned on the
  host. Fix: `--user (id -u):(id -g)`, paired with `-e HOME=/tmp` since
  the host UID doesn't exist in the container's `/etc/passwd`.

- **2026-09-04** — Current case parameters give **Re = U·L/ν =
  1×1/1e-5 = 100,000** (`constant/physicalProperties`: `nu 1e-05`;
  `0/U`: lid speed 1 m/s; cavity side length 1 m) — far above the
  ~10,000 upper bound for laminar lid-driven-cavity benchmarks (Ghia et
  al. 1982). `constant/momentumTransport` has `simulationType laminar`
  with a `RAS { model kEpsilon; turbulence on; }` block underneath it
  that is entirely dead code while `simulationType` stays `laminar`.
  Flagged to the user; **left unchanged** — user is adjusting Re
  themselves.

- **2026-09-04** — `system/controlDict`'s fixed `deltaT 0.05` (no
  `adjustTimeStep`) combined with the mesh's cell size (~1/49 ≈ 0.02 m,
  50-point transfinite edges) and lid speed 1 m/s gives Courant number
  ≈ 2.4 — above the usual CFL ≈ 1 guideline, especially with `PIMPLE`'s
  `nCorrectors 2` and no `nOuterCorrectors` set (defaults to 1, i.e.
  effectively PISO-like, less tolerant of large Courant numbers than a
  true outer-corrected PIMPLE loop). Flagged; **left unchanged** — user
  is adjusting the time step themselves.

- **2026-09-04** — `system/fvSolution`'s `PIMPLE` block already sets
  `pRefCell 0; pRefValue 0;`, which resolves the classic all-Neumann
  pressure singularity for a fully closed cavity (no outlet/inlet, every
  boundary is `zeroGradient` or `empty`). Worth checking for any closed-
  domain incompressible case — without a pressure reference, the solve
  is only defined up to an arbitrary constant.

- **2026-09-07** — OpenFOAM's own `sample` function object
  (`system/sampleDict`, `type sets;`, run via `postProcess -func
  sampleDict`) does cell-to-point interpolation along an arbitrary line
  using the mesh's own connectivity — the right tool for extracting a
  `y=const`/`x=const` profile, rather than trying to reconstruct a
  spatial grid from `internalField`'s raw per-cell row order (which
  follows OpenFOAM's internal cell numbering from `gmshToFoam`, not a
  simple i,j scan — confirmed by inspecting `1/U`: `nonuniform
  List<vector>` with 2401 entries, i.e. per-*cell* not per-node, in
  mesh-numbering order).

- **2026-09-07** — Reference-vs-present-result plotting convention for a
  validation figure: tabulated benchmark data (few discrete points) as
  markers only, the current simulation (densely-sampled) as a line only.
  Why it matters: the original `plots.gp` drew both `with line`, which
  visually implies Ghia's 17 discrete grid values are a continuous
  measurement — misleading for what's actually a coarse reference table.

- **2026-09-07** — `git status`'s rename detection (`R` entries) pairs a
  staged deletion with a staged addition purely by *content similarity*,
  not by any actual `mv`/`git mv` relationship — e.g. it paired a deleted
  `cylinderLowRe/constant/physicalProperties` with an unrelated, newly
  added `cavity/openFoam/1/uniform/time` file just because both happened
  to be short/similar text. Why it matters: don't trust an `R` line in
  `git status` as proof a file was intentionally moved — cross-check with
  `find`/`git log` on the actual path before concluding a "move" was a
  move and not a coincidental deletion+addition pair.

- **2026-09-07** — Untracked files don't travel with a git-aware
  directory reorg: `plot_notes.md` and `plots.gp.bak` existed on disk but
  weren't yet `git add`ed when the user's reorg ran, so they were left
  behind at the old path while every *tracked* file in the same directory
  moved cleanly (showing as clean `R` renames in `git status`). Why it
  matters: after any large restructuring, check for untracked files at
  the old location specifically — they're invisible to a `git mv`-based
  move and to `git status`'s rename pairing, so they don't show up as
  "something went wrong" the way a tracked file would.

- **2026-09-07** — A `.gitignore` pattern anchored to a literal path
  prefix (e.g. `navStokes/openFoam/**/...`) fails silently on a directory
  rename — no error, no warning, the previously-ignored files just quietly
  become trackable again. Confirmed via `git check-ignore -v` before/after
  the fix. Anchoring on a path *segment name* instead (`navStokes/**/
  openFoam/**/...`, matching `openFoam` at any depth) is robust to the
  container/subdirectory relationship flipping, as it did here.

## Future Action Items

- [ ] `Allrun` still runs `blockMesh` + `mirrorMesh` (copied from the
      `cylinderLowRe`/tutorial template — that case itself has since been
      deleted, 2026-09-07) rather than the actual `mesh.py` (gmsh) +
      `gmshToFoam` workflow this case now uses — noticed while reviewing
      the case this session, not yet fixed.
- [x] User to pick a target Reynolds number (via `constant/
      physicalProperties` `nu` and/or lid speed in `0/U`) and set
      `system/controlDict`'s `deltaT`/`adjustTimeStep` accordingly — done
      2026-09-04, see Status above.
- [ ] Consider whether the `lid` physical surface should be retyped so
      `gmshToFoam` assigns it `type wall` (matching physical reality)
      instead of the current generic `type patch`, if wall-function/
      force/`yPlus` post-processing is ever needed on it.
- [ ] Run currently stops at `t=161` of `endTime 200` — not confirmed in
      this session whether that's a completed/converged stopping point or
      the run was just paused; worth checking before treating the
      `postProcessing/sampleDict/161/` profiles as final.
- [ ] `plots.gp`'s output is still raster PNG; vector output (`pdfcairo`/
      `epslatex`) considered for a journal-ready figure but not applied —
      see "Options considered but not applied" in `plot_notes.md`.
- [ ] Today's `.gitignore` fix and the reorg itself are still uncommitted
      (working tree only) — review `git status`/`git diff` before
      committing, given how much moved in this session.

## References

- Ghia, U., Ghia, K.N., Shin, C.T. (1982), *High-Re solutions for
  incompressible flow using the Navier-Stokes equations and a multigrid
  method*, J. Comput. Phys. 48(3):387–411. (Benchmark laminar
  lid-driven-cavity solutions, Re up to 10,000 — basis for the Re
  ceiling noted above.)
