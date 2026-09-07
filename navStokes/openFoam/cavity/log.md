# navStokes/openFoam/cavity — Log

_Started: 2026-09-04 | Last updated: 2026-09-04_

## Status

- **2026-09-04** — Re and `deltaT` findings below (flagged, left
  unchanged at the time) have since been fixed by the user: `nu` is now
  `1e-02` (Re=100, standard Ghia et al. laminar benchmark case) and
  `deltaT` is `0.01` (Courant ≈0.49, comfortably under CFL≈1). Re-checked
  and confirmed — see Learnings for the numbers this superseded.

## Objective

OpenFOAM case for classic 2D lid-driven cavity flow — a unit square, top
("lid") wall moving tangentially at a fixed speed, other three walls
stationary. Sibling case to `../channel_flow/` under `navStokes/openFoam/`.
Mesh is built with a custom `mesh.py` (gmsh Python API) rather than
`blockMesh` or a `.geo` script, converted to OpenFOAM's `constant/polyMesh`
via `gmshToFoam`, and run with `foamRun -solver incompressibleFluid`.

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

## Future Action Items

- [ ] `Allrun` still runs `blockMesh` + `mirrorMesh` (copied from the
      `cylinderLowRe`/tutorial template) rather than the actual
      `mesh.py` (gmsh) + `gmshToFoam` workflow this case now uses —
      noticed while reviewing the case this session, not yet fixed.
- [ ] User to pick a target Reynolds number (via `constant/
      physicalProperties` `nu` and/or lid speed in `0/U`) and set
      `system/controlDict`'s `deltaT`/`adjustTimeStep` accordingly (see
      Learnings above for the numbers this session found).
- [ ] Consider whether the `lid` physical surface should be retyped so
      `gmshToFoam` assigns it `type wall` (matching physical reality)
      instead of the current generic `type patch`, if wall-function/
      force/`yPlus` post-processing is ever needed on it.

## References

- Ghia, U., Ghia, K.N., Shin, C.T. (1982), *High-Re solutions for
  incompressible flow using the Navier-Stokes equations and a multigrid
  method*, J. Comput. Phys. 48(3):387–411. (Benchmark laminar
  lid-driven-cavity solutions, Re up to 10,000 — basis for the Re
  ceiling noted above.)
