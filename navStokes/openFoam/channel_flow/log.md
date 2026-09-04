# navStokes/openFoam — Log

_Started: 2026-09-02 | Last updated: 2026-09-03_

## Objective

OpenFOAM baseline case for a 3D turbulent channel flow between two infinite
parallel plates, to validate the FEM (dolfinx) Navier-Stokes solver being
built in `navStokes/` against. See `../overview.md` for the full project
framing (objectives, assumptions, FEM build order).

## Decisions

- **2026-09-02** — OpenFOAM baseline uses the openfoam.org (Foundation)
  line, run via the `microfluidica/openfoam:org` Docker image, which is
  actually **v14** (checked `$WM_PROJECT_VERSION` inside the container,
  not v13 as originally assumed).
  Why: user has an existing Docker setup for this line; v14 replaces the
  originally-planned v13.
  Ruled out: openfoam.com (ESI/OpenCFD) line — different solver/case
  conventions, would fragment the comparison.

- **2026-09-02** — Case dictionaries use the modular `physicalProperties` /
  `momentumTransport` naming (the `foamRun -solver incompressibleFluid`
  convention), not the legacy `transportProperties` / `turbulenceProperties`.
  Why: confirmed via the container's own `incompressibleFluid/pitzDaily`
  tutorial — v11+ of the Foundation line modularized solvers into
  `foamRun`; standalone `icoFoam`/`pimpleFoam` apps are not the current
  convention.

- **2026-09-02** — `.gitignore` OpenFOAM rules (mesh, time-step dirs,
  `processor*/`, `postProcessing/`, `VTK/`, logs) are scoped to
  `navStokes/openFoam/**` rather than repo-wide.
  Why: avoid accidentally ignoring unrelated files elsewhere in the repo
  (e.g. any future numeric-named directories in the 1D/2D FEM projects).

- **2026-09-02** — Channel geometry spans the full height (walls at
  y=-1 and y=+1), not a half-channel with a symmetry plane at y=0.
  Why: simpler physical-group set, no reliance on a symmetry BC assumption
  for what's meant to be a genuinely turbulent (asymmetric-instant, only
  statistically symmetric) flow.
  Ruled out: half-channel + symmetry-at-centerline, which was the first
  draft (physical group `pipe_symmetry` / `centerline`) but got dropped.

- **2026-09-02** — Volume mesh is structured hexahedra (`Transfinite
  Curve`/`Surface` + `Recombine` + `Extrude{...Layers{N}; Recombine;}`),
  not an unstructured tet mesh.
  Why: practical DNS/LES-scale channel resolution needs a hex mesh with
  wall-normal grading; an unstructured tet mesh at the resolution DNS
  requires would be on the order of 10^8 elements (checked before
  committing to this approach).

- **2026-09-03** — Target Re_τ = 180, matching Kim, Moin & Moser (1987)
  directly (Re_centerline = 3300, ν computed for air at ρ=1.225 kg/m³,
  μ=1.789e-05 Pa·s). First-cell **center** placed at y⁺ = 0.05, per
  KMM's own reported near-wall resolution.
  Why: Re_τ = 180 is KMM's *measured* DNS result for this Re, not an
  estimate — more trustworthy than deriving u_τ from a generic
  correlation when the point is to match their benchmark exactly.
  Ruled out: estimating u_τ via the Schlichting flat-plate skin-friction
  formula, Cf=(2·log10(Re)−0.65)^-2.3, applied with L=h and U=U_c. This
  is an *external* boundary-layer correlation, not calibrated for
  internal/duct flow — it implied Re_τ≈277 here, ~54% off KMM's actual
  180.

- **2026-09-03** — Case is run as `simulationType laminar;` in
  `constant/momentumTransport`, not RAS/LES, with no turbulence-model
  fields in `0/` (`k`, `epsilon`, `omega`, `nut`, `nuTilda` removed).
  Why: this is a DNS run — the point is to resolve the Navier-Stokes
  equations directly at the mesh's own resolution, so no turbulence
  model or wall function is used; the case is set up exactly as if the
  flow were laminar, and turbulence emerges from the resolved solution.

## Tasks

- [x] `../overview.md` written — objectives, assumptions, decisions, FEM
      build order, open risks (2026-09-02)
- [x] `.gitignore` updated for OpenFOAM outputs, scoped to
      `navStokes/openFoam/` (2026-09-02)
- [x] `mesh.geo` — geometry, curve loop, and physical groups (`wall`,
      `left`, `right`, `inner`, `outer`, `domain`) fixed and verified
      (syntax bugs, stale surface-tag references, curve-loop
      renumbering) (2026-09-02)
- [x] `mesh.geo` — structured hex volume mesh via `Layers{}`/`Recombine`
      in the `Extrude`, verified 0 tets/pyramids (2026-09-02)
- [x] `mesh.geo` — wall-normal `Bump` grading direction bug found and
      fixed (2026-09-02)
- [ ] `createOpenFoam.fish` — still incomplete. Currently only creates
      `system/{controlDict,fvSchemes,fvSolution}` (empty) and
      `constant/polyMesh/*` (empty placeholders) plus a stray, incorrect
      `constant/Properties/` directory. Still needs: `constant/
      physicalProperties`, `constant/momentumTransport`, `0/U`, `0/p`
      (FoamFile header only, no parameters, so it's reusable across
      future projects) — this was interrupted mid-task and never
      finished.
- [x] Pick a target Re_τ for the channel case and re-derive `mesh.geo`'s
      wall-normal resolution from it (2026-09-03) — Re_τ=180 (KMM 1987),
      first-cell-center at y⁺=0.05.
- [x] Re-derive `mesh.geo`'s streamwise/spanwise resolution (Δx⁺, Δz⁺)
      for Re_τ=180 (2026-09-03) — grid is now 192×129×160 points
      (191×128×159 divisions): Δx⁺=11.84 (target ~10-15), Δz⁺=7.11
      (target ~5-7, right at the edge), first-cell-center y⁺=0.0500
      exactly. ~4.0M nodes / ~3.66M hex elements.
- [x] Convert `mesh.msh` to OpenFOAM's `constant/polyMesh` via
      `gmshToFoam` (2026-09-03) — required re-exporting the mesh as MSH
      format 2.2 first (see Learnings); 3,887,232 hex cells, all 5
      patches mapped correctly.
- [x] Case set up as `simulationType laminar` (DNS — no turbulence
      model/wall functions), turbulence-model fields removed from `0/`
      (2026-09-03).
- [ ] Start the FEM build order from `../overview.md` (steady 2D Stokes +
      strong BCs + Taylor-Hood, verified against the analytic Poiseuille
      profile) — not yet started.

## Learnings

- **2026-09-02** — Gmsh's `Extrude {dx,dy,dz} { Surface{tag}; };` returns
  `out[0]` = far cap, `out[1]` = volume, `out[2..]` = side surfaces in
  `Curve Loop` order.
  Why it matters: don't trust this mapping from memory or from a comment
  in the file — verify per-geometry with `gmsh.model.getBoundingBox()`.
  Renumbering points/lines changes the mapping, and a stale comment
  describing the old mapping looks identical to a correct one.

- **2026-09-02** — `Transfinite Curve{...} Using Bump C`: `C > 1` clusters
  nodes at the curve's **center**; `C < 1` clusters at the **endpoints**.
  Why it matters: this is the opposite of the intuitive reading of "bump."
  For wall clustering in a channel mesh, always use `C < 1`. Verified
  empirically on an isolated single-line test before trusting it on the
  real geometry — caught a real bug this way (`Bump 1.2`/`1.4` had been
  clustering at the centerline, not the walls).

- **2026-09-02** — Gmsh's CLI dispatches how it parses an input file
  purely by its file extension: `.msh` is read as mesh *data*, not a
  script. A file containing `.geo`-script syntax but literally named
  `mesh.msh` fails to parse (or silently produces an empty mesh) under
  `gmsh -3 mesh.msh`.
  Why it matters: geometry scripts belong in `.geo` files; the generated
  mesh output belongs in `.msh`. (This is why the file is now `mesh.geo`,
  self-saving to `mesh.msh` via an embedded `Save "mesh.msh";`.)

- **2026-09-02** — `Recombine Surface{tag};` only quad-recombines the 2D
  base surface. The 3D `Extrude{...}` needs its own `Layers{N};
  Recombine;` *inside* the extrude block to produce hexahedra — without
  it, extruding a quad surface still fills the volume with a hybrid
  tet/pyramid mesh.
  Why it matters: this is the actual fix for "a proper hex mesh," not
  just recombining the base surface (confirmed by checking element type
  counts before/after: hybrid tet/pyramid mesh had 0 hexahedra type-5
  elements; after adding `Layers{}`/`Recombine{}` to the extrude, it was
  100% type-5 hexahedra + type-3 quad boundary faces).

- **2026-09-02** — DNS near-wall resolution target: first-cell y⁺ ≲ 1,
  Δx⁺ ≈ 10-15 (streamwise), Δz⁺ ≈ 5-7 (spanwise) — Lee & Moser (2015);
  Kim, Moin & Moser (1987). RANS: y⁺ ~ 1 for low-Re/wall-resolved models,
  y⁺ ∈ [30, 300] for wall functions (avoid the 5-30 buffer-layer gap).
  Why it matters: Δy_wall = y⁺_target · h / Re_τ converts a target y⁺
  into an actual mesh spacing once Re_τ is picked. Checked the *current*
  `mesh.geo` against this: first-cell spacing is 0.0177 (h=1), which
  implies y⁺₁ ≈ 3.2 at Re_τ=180 — too coarse for DNS, borderline even
  for wall-resolved RANS. Streamwise/spanwise resolution is currently
  even coarser relative to the DNS target at that Re_τ.

- **2026-09-03** — y⁺=0.05 at Re_τ=180, h=1 gives a wall distance of
  y = y⁺·h/Re_τ = 2.778e-4 m (this simplifies to y⁺·h/Re_τ directly —
  ν cancels once u_τ is written as Re_τ·ν/h, so ν itself doesn't need
  to appear in the final formula).
  Why it matters: caught two arithmetic slips working toward this
  number — using the Schlichting flat-plate correlation instead of
  KMM's actual Re_τ (see Decisions above), and separately using μ where
  ν was needed (u_τ = Re_τ·ν/h, not Re_τ·μ/h) — before landing on this
  value. Since this is the *cell-center* target, not the first-node
  spacing, the actual node placed off the wall needs to be at 2×y =
  5.556e-4 m (the cell center sits halfway between the wall face and
  the first interior node).

- **2026-09-03** — Retuned `Transfinite Curve{12,14}` (wall-normal) from
  40 pts/Bump 0.2 to **120 pts/Bump 0.01054**, found by bisecting the
  Bump coefficient on an isolated single-line test (`bump_test.geo` in
  scratch) until the first node landed at 5.556e-4 m. Verified on the
  real geometry: first spacing 5.5543e-4 m (target 5.556e-4), growth
  ratio ~1.10 near the wall, mesh still 100% hex (92820 type-5 elements,
  15602 type-3 boundary quads, 0 tets/pyramids).
  Why it matters: growth ratio matters as much as the first-cell size —
  tried 80/100/120/150/200 pts for the same target spacing and growth
  ranged 1.18 (80 pts) down to 1.05 (200 pts); picked 120 (~1.10) as a
  balance between mesh size and cell-to-cell size jump near the wall.
  Streamwise and spanwise resolution were redone the same day — see the
  next two entries; superseded, don't use the N=120/Bump 0.01054 numbers
  above as current.

- **2026-09-03** — `Transfinite Curve{c} = N` sets N *points* along the
  curve (N-1 divisions/cells) — confirmed empirically (`Transfinite
  Curve{1}=192` on an isolated line produced 192 nodes, 191 line
  elements). `Extrude {...} { Surface{s}; Layers{N}; Recombine; }` is
  different: `Layers{N}` sets N *divisions* directly (N+1 points) on the
  extruded direction, with no separate `Transfinite Curve` needed or
  wanted on the spanwise edges it creates.
  Why it matters: `Transfinite Curve{11,13} = 192` was silently giving
  191 x-divisions, not 192 — an off-by-one that's easy to carry through
  a whole mesh-sizing calc undetected. The two commands (`Transfinite
  Curve` vs. extrude `Layers{}`) use different counting conventions for
  the same-looking integer argument; check which one before trusting an
  intended division count.

- **2026-09-03** — Final grid: 192×129×160 points (191×128×159
  divisions) → Δx⁺=11.84, Δz⁺=7.11, first-cell-center y⁺=0.0500, for
  Re_τ=180. `Transfinite Curve{12,14}` (wall-normal) is 129 pts/Bump
  0.011566 (retuned from the 120-pt value above — the Bump coefficient
  needed for a given first-cell spacing depends on point count, so it
  has to be re-bisected whenever N changes, not just carried over).
  ~4.0M nodes / ~3.66M hex elements.

- **2026-09-03** — To check min/max element edge length on a *structured,
  axis-aligned* mesh, no gmsh Python API is needed (and the `gmsh` pip
  module wasn't installed here — only the CLI binary): just collect the
  distinct node coordinates along each axis from the `.msh` file, take
  consecutive differences, and min/max those. Every element edge on this
  mesh is exactly a dx, dy, or dz spacing since the grid is axis-aligned.
  Result: min edge 5.5561e-4 m (wall-normal, at the wall), max edge
  6.5759e-2 m (uniform in x), ratio ~118 — a large but expected spread
  for a wall-resolved mesh (near-wall cells are short in y, long in x/z,
  matching the elongated near-wall turbulent structures).
  Why it matters: for a genuinely unstructured or non-axis-aligned mesh
  this shortcut doesn't work — would need actual element connectivity
  (via `pip install gmsh` and `gmsh.model.mesh.getElements`) to walk
  each element's edges individually.

- **2026-09-03** — `gmshToFoam` (openfoam.org / Foundation line, v14)
  fails on gmsh's default output format. gmsh 4.15.2 writes MSH format
  4.1 by default; `gmshToFoam` mis-parses its `$Entities`/`$Nodes`
  block structure and crashes (`FOAM FATAL IO ERROR: Attempt to get
  back from bad stream`) after misreading the entity count (27, i.e.
  8 points + 12 curves + 6 surfaces + 1 volume) as the node count.
  Reproduced directly in the `microfluidica/openfoam:org` container
  before touching anything, to confirm root cause rather than guess.
  Fix: export MSH format **2.2** instead — `gmsh ... -format msh2` (now
  set permanently via `Mesh.MshFileVersion = 2.2;` in `mesh.geo`, right
  before `Save`). Re-ran `gmshToFoam` on the 2.2 output: clean read,
  3,887,232 hex cells (matches 191×128×159 divisions exactly), all 5
  patches (`inner`,`wall`,`right`,`left`,`outer`) mapped correctly.
  Why it matters: **for any future project using openfoam.org's
  `gmshToFoam`, always save gmsh meshes as MSH format 2.2, not gmsh's
  current default (4.1).** This isn't specific to this mesh or
  geometry — it's a `gmshToFoam`/gmsh-version compatibility gap that
  will resurface on any new case built the same way.

## Future Action Items

- [ ] Finish `createOpenFoam.fish` / the minimal empty-dictionary case
      skeleton (`physicalProperties`, `momentumTransport`, `0/U`, `0/p`,
      header-only so it's reusable for future projects). Also check
      whether openfoam.org ships any built-in "new case" command before
      hand-rolling further — as of v14, `foamCloneCase` exists but only
      clones from an existing source case/template; there is no built-in
      truly-empty-case generator.
- [ ] Decide the `gmshToFoam` (or alternative) mesh-conversion path, and
      set up `inner`/`outer` as a `cyclic` boundary-condition pair for
      the periodic spanwise direction.
- [ ] Begin the FEM build order from `../overview.md`: steady 2D Stokes +
      strong BCs + Taylor-Hood, verified against the analytic Poiseuille
      profile, before attempting the full 3D turbulent Nitsche solver.

## References

- Kim, Moin & Moser (1987), *Turbulence statistics in fully developed
  channel flow at low Reynolds number*, JFM 177:133–166.
- Moser, Kim & Mansour (1999), *DNS of turbulent channel flow up to
  Re_τ=590*, Phys. Fluids 11:943–945.
- Lee & Moser (2015), *Direct numerical simulation of turbulent channel
  flow up to Re_τ≈5200*, JFM 774:395–415 ([arXiv:1410.7809](https://arxiv.org/abs/1410.7809)).
- Jiménez & Moin (1991), *The minimal flow unit in near-wall turbulence*,
  JFM 225:213–240.
- Wilcox, D.C., *Turbulence Modeling for CFD* (3rd ed., 2006), DCW
  Industries.
- Pope, S.B., *Turbulent Flows* (2000), Cambridge University Press.
- Menter, F.R. (1994), *Two-Equation Eddy-Viscosity Turbulence Models for
  Engineering Applications*, AIAA Journal 32(8):1598–1605.
- [Oden Institute channel-flow DNS database](https://turbulence.oden.utexas.edu/channel2015/content/README_2015.html)
  (Re_τ = 180, 550, 1000, 2000, 5200 — usable as a second validation
  baseline alongside the OpenFOAM run).
