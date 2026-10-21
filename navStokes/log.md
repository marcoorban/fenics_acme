# navStokes — Project Log

_Started: 2026-09-09_

This is the **project-wide** log: tooling, workflow, and technical gotchas that
apply across cases, or that were learned in one case but generalize to any future
one. Case-specific decisions/tasks (mesh parameters, boundary conditions, Reynolds
numbers, benchmark comparisons for that specific geometry) stay in the per-case logs:

- `cavity/log.md` — lid-driven cavity case.
- `cylinder/log.md` — flow past a cylinder (2D DFG 2D-1/2-2/2-3 benchmarks complete; 3D,
  square and circular cross-sections, is next — see the new `#openFoam` entry below).
- `cavity/plot_notes.md` — detailed rationale for the `plots.gp` styling pass (kept
  separate since it's long and cavity-specific in its examples, even though the
  *convention* it establishes is project-wide — synthesized into `#plotting` below).

## How entries here are tagged

Each entry has a **scope** tag and one or more **topic** tags:

- **Scope** — `[project]` (repo/tooling-wide, not about any one case) or
  `[cavity]`/`[cylinder]` (the lesson was learned working on that case, but
  is written up here because it's expected to generalize to other cases too).
- **Topic** — `#openFoam`, `#gmsh`, `#git`, `#docker`,
  `#plotting`, `#fem` — the subsystem the note is actually about.

## #openFoam

- **[project]** — **2026-09-14: 2D cylinder flow simulations are done.**
  All three DFG benchmark cases are complete and match published reference
  values closely — `2D-1` (steady, `Re=20`), `2D-2` (unsteady vortex
  shedding, `Re=100`), `2D-3` (time-ramped inlet, `Re=100`, `sin(pi*t/8)`
  forcing). Full numbers, methods, and plots are in `cylinder/log.md` and
  `cylinder/results.csv` — nothing case-specific repeated here. **Next step: 3D
  cylinder flow simulations, with both square and circular cross-sections.**

- **[cavity]** — `internalField` in a field file (`U`, `p`, ...) is per-*cell*,
  in the mesh's own internal numbering — not a reshape-able `i,j` grid, even
  on a structured hex mesh. For a profile along `x=const`/`y=const`, use
  the `sample` function object (`system/sampleDict`, `type sets;`, run via
  `postProcess -func sampleDict`) — it interpolates via the mesh's actual
  connectivity (`interpolationScheme cellPoint`) rather than assuming a scan order.

- **[cylinder]** — gmsh's own `"N nodes M elements"` summary bundles 2D surface/boundary
  elements with 3D cells — not a usable proxy for OpenFOAM's actual cell count,
  which can be anywhere from ~0.33x to ~0.95x of that total depending on resolution.
  Check `checkMesh`'s `cells:` line instead when targeting a specific cell
  count, and expect to iterate empirically rather than trust `1/size²` scaling
  once a background sizing field with a fixed transition distance is involved.

- **[cylinder]** — For a pseudo-2D single-layer extrusion, always add
  `recombine=True` to `gmsh.model.geo.extrude(...)` (as `cavity`'s
  script always did) — without it, gmsh tetrahedralizes/subdivides the
  layer instead of producing clean prisms, which measurably worsens mesh
  non-orthogonality/skewness for no benefit on a geometry this simple.

- **[cylinder]** — `foamPostProcess -func <name>` needs an explicit `-solver
  <name>` flag when the function object (e.g. `forces`/`forceCoeffs`) looks
  up solver fields (`U`, `p`) — without it, the mesh and time directory
  load fine but the fields are never registered, giving a confusing
  `Could not find U, p` even when the files genuinely exist on disk.

- **[cylinder]** — `gmshToFoam`'s automatic `wall`/`empty` patch-type detection
  (documented below as working via physical-surface name) is **not consistently
  reliable** — it worked on one mesh revision and silently fell back to generic
  `type patch;` for every patch (including `empty`, which is fatal downstream) on
  another, with the exact same physical-group names and no warning at conversion
  time. **Don't trust it — verify `constant/polyMesh/boundary` after every mesh
  regeneration**, and force the types explicitly rather than relying on the
  heuristic: `changeDictionary` is deprecated (removed from this project 2026-09-14)
  — use `foamDictionary -set` on `constant/polyMesh/boundary` instead (see
  `cylinder/openFoam/mesh.sh` for the pattern: `gmshToFoam` then `foamDictionary -set`).

- **[cylinder]** — Before diagnosing a process as "hung" (e.g. an MPI deadlock
  after a partial crash), check actual process state (`docker top`/`ps`) rather
  than just a long container uptime — `docker ps` showing a container "Up 3 days"
  turned out to be an unrelated idle interactive shell, not the run in question
  (which had already exited and auto-removed itself via `docker run --rm`).

- **🚨 [cylinder]** — When retesting a fix, change **one variable at a time**!! After
  fixing an `nNonOrthogonalCorrectors`-related blow-up, the retest also jumped
  `deltaT` 100× in the same run — so when it hit a *different* Courant-number failure
  (mean elevated too, not just max — a plain CFL issue this time, see `[cylinder]`
  log for the full mean-vs-max diagnostic), there was no way to tell whether the
  original fix actually worked or not. Bundling "does my fix work" with "how far
  can I push this other setting" into one test destroys attribution either way.

- **🚨 [cylinder]** — A "Courant number blowing up" crash is **not automatically a
  timestep problem**, even though it looks exactly like one!! Before shrinking `deltaT`,
  check **mean vs. max** Courant number separately in the log — they diagnose two
  different failure modes: mean too high domain-wide → genuine CFL/timestep issue,
  smaller `deltaT` actually helps; mean pinned near-zero while max explodes at one
  localized spot → numerical instability (bad cell / insufficient non-orthogonal
  correction), and shrinking `deltaT` does **not** fix it — confirmed directly:
  dropping `deltaT` `1e-5 → 1e-6` still blew up at the identical physical time. In
  the end the fix that actually worked wasn't the correctors — it was regrading the
  mesh itself (a `Distance`+`Threshold` field for smooth sizing, `recombine=True` for
  clean prisms — see `cylinder/log.md` for the full writeup); a coarse pass on the
  regraded mesh was confirmed stable before refining cell count for the real run.

- **[cylinder]** — `Allclean` in both cases only removed top-level time directories
  — it never touched `processorN/`!! Those directories carry their *own* time
  directories and their *own* mesh addressing (`cellProcAddressing` etc.), both
  sized for whatever mesh existed when `decomposePar` last ran. Regenerate the mesh
  afterward (e.g. a `mesh.py` parameter change) and those leftovers don't error at
  `decomposePar` time — they sit there silently until `reconstructPar -latestTime` (or
  a solve reaching a coincidentally matching write time) picks up the stale data and
  collides with the freshly-rebuilt addressing: `FOAM FATAL IO ERROR: size X is not
  equal to the given value of Y`!! Caught in `cylinder` after several mesh revisions
  left `processor0-7` holding a full old solved run (~2080 cells/subdomain) while
  the current mesh decomposes to ~13352/subdomain. Fixed **both** `Allclean` scripts
  to also `rm -rf processor*/` — same latent bug was sitting in `cavity` too, just
  hadn't been triggered yet since its mesh hasn't changed since its last parallel run.

- **[cavity]** — `gmshToFoam` infers a boundary's OpenFOAM `type` (`wall` vs.
  generic `patch`) heuristically from the *physical-surface name* — a group
  literally named `"wall"` gets `type wall`; one named `"lid"` (despite being a
  real, moving wall) gets `type patch`. Harmless for a laminar Dirichlet velocity
  BC, but would matter for wall-function/`wallShearStress`/`yPlus` post-processing
  on that patch. Worth checking on any future case's `constant/polyMesh/boundary`.

- **[project]** — `fvSolution`'s `solvers{}` `tolerance`/`relTol`
  and `PIMPLE{}`/`SIMPLE{}`'s `residualControl` are **not** the
  same kind of "convergence" setting, easy as that is to assume:
  - `tolerance`/`relTol` (per solver, e.g. `p`, `U`) control how tightly that *one
    linear system* is solved **within a single time step** — numerical accuracy of that
    step's solve, nothing to do with whether the physical solution has stopped evolving.
  - `residualControl` (in `PIMPLE`/`SIMPLE`) is the actual steady-state detector:
    it checks each listed field's *initial* residual every time step and auto-stops
    the run once all of them are below threshold simultaneously — the right
    tool for "stop once it's converged" instead of always running to `endTime`.

## #plotting

Full rationale and options considered live in `cavity/plot_notes.md`;
synthesized here as the **project's plotting convention** (applies to any
future case's validation plots, e.g. cylinder drag/lift/Strouhal comparisons):

- **[cavity]** — Before plotting a *raw* OpenFOAM field value against a *published
  benchmark* table, check whether the benchmark is non-dimensionalized. Ghia et al.
  (1982)'s `u`/`v` are scaled by the lid speed (`u* = u/U_lid`) — invisible when
  `U_lid=1` (numerically identical to the raw value), but a stark mismatch once
  `U_lid` changes (e.g. Re=400 run uses `U_lid=4`: raw `U_x` ranges to 4.0, Ghia's
  column ranges to 1.0 — same shape, 4× the scale). Verified the normalization
  itself directly from the data: the lid-boundary row (`y=1`) is exactly `1.00000`
  across *every* Reynolds-number column in Ghia's table, which only holds if
  velocities are reported as a fraction of the (fixed, non-dimensional) lid speed.
  General lesson: check a benchmark table's boundary-condition row for this kind
  of tell before assuming raw and tabulated values are directly comparable.

- **[project]** — `pngcairo enhanced font "Helvetica,14"` + italicized axis
  labels (`{/Italic y}`) — matches ISO 80000 / journal-style convention for
  physical-quantity symbols. Vector output (`pdfcairo`/`epslatex`) considered for
  a genuinely journal-ready figure but not applied yet — raster PNG kept for now.

- **[project]** — Okabe-Ito colorblind-safe palette (black markers +
  `#0072B2` blue line) instead of gnuplot's default red/green cycling.

- **[project]** — Validation-figure convention: tabulated *benchmark*
  data (few discrete points) as markers only; the *current simulation*
  (densely sampled) as a line only. Both `with line` misleadingly
  implies a coarse reference table is a continuous measurement.

## #gmsh

- **[cavity]** — `gmsh.model.geo.extrude(...)` (Python API) returns entities
  as `[far cap, volume, lateral surfaces in Curve Loop order]` — same
  convention as the `.geo`-script `Extrude{}` command. Verify per-geometry
  against the resulting `constant/polyMesh/boundary` face counts rather than
  trusting the mapping from memory — renumbering points/lines changes it.

## #git

- **[project]** — A delete/modify merge conflict (`deleted by us`) means your
  branch removed a file while the other branch kept/changed it — git checks "their"
  version into the working tree with no conflict markers (nothing to textually
  merge), but still requires an explicit `git add`/`git rm` to resolve. If the
  deletion was intentional, `git rm <path>` confirms it and stages the resolution.

- **[project]** — Untracked files don't travel with a git-aware directory reorg
  — they're invisible to `git mv`/rename-pairing, so they silently get left
  behind at the old path while every tracked sibling moves cleanly. Check for
  untracked files at the old location specifically after any large restructuring.

- **[project]** — `git status`'s `R` (rename) lines are a *content
  similarity* heuristic, not proof of an actual `mv`/`git mv`. Seen
  pairing an unrelated staged deletion with an unrelated staged addition
  purely because both files were short and similar. Cross-check with
  `find`/`git log` on the real paths before trusting an `R` line.

- **[project]** — `.gitignore`'s OpenFOAM section is anchored on the path
  *segment name* `openFoam` at any depth (`navStokes/**/openFoam/**/...`), not a
  literal prefix (`navStokes/openFoam/**/...`). A literal-prefix anchor breaks
  silently (no error, previously-ignored files just become trackable again) the
  moment a case's directory layout changes — which happened once already (case
  moved from `navStokes/openFoam/<case>/` to `navStokes/<case>/openFoam/`,
  ~653 generated files went from ignored to staged before the fix).

## #docker

- **[project]** — "permission denied" talking to `/var/run/docker.sock` after adding
  yourself to the `docker` group usually means the *current shell session* hasn't
  picked up the new group membership yet (group lists are fixed at login/process
  creation) — confirm with `id`/`groups` showing `docker` missing, then fix with `newgrp
  docker` (immediate, current shell only) or a fresh login/terminal (persists). Run
  in your own terminal, not through a tool that spawns isolated one-off shells per
  command — the group change won't carry over between separate invocations there.

- **[project]** — `oforg.fish` (the Docker wrapper for `microfluidica/openfoam:org`)
  must run as `--user (id -u):(id -g) -e HOME=/tmp`, not root. `-v $PWD:/project`
  bind mounts share the host's user namespace, so a container run without
  `--user` writes root-owned output back to the host; `HOME=/tmp` is needed
  alongside it because the host UID has no entry in the container's `/etc/passwd`.

## #fem

- **[cavity]** — **The Navier-Stokes convective term `(u.grad)u` is nonlinear
  because the unknown plays two roles at once — it's both the transported
  quantity and the transport velocity.** Contrast with a passive
  advection-diffusion problem (`-D*grad^2(f) + a.grad(f) = source`): there,
  the transport velocity `a` is externally given data, so the equation is
  linear in the one unknown `f` and assembles into a single `LinearProblem`
  with no iteration needed. In momentum transport, `a` *is* `u` — the thing
  being solved for — so `(u.grad)u` is quadratic in the unknown and can't be
  assembled into a bilinear form directly.
  **Fix: Picard linearization** — replace *one* copy of `u` in the
  convective term with a known/lagged field `u_n` (playing `a`'s role),
  and keep the *other* copy as the real trial function `u` (playing `f`'s
  role): `inner(dot(u_n, nabla_grad(u)), w)*dx`. This restores the passive-
  scalar structure, so the term is linear in `u` for a *given* `u_n` and
  belongs in the system matrix, not the right-hand side. Iterate — solve,
  set `u_n <- u`, resolve — until `u_n` and `u` agree to within tolerance;
  at that fixed point the Picard-linearized equation coincides with the
  true nonlinear steady NS equation.
  **What NOT to do (caught in `cavity/fem/fem_strong.py` and the two
  `fem_weak*.py` variants, fixed 2026-10-21):** the original scheme froze
  *both* copies of `u` at `u_n` (`inner(dot(u_n, nabla_grad(u_n)), w)*dx`)
  and moved the entire term to the RHS as a precomputed forcing — the
  trial function `u` never appeared in the convective term at all. This
  is NOT Picard; it's a single fully-explicit/Jacobi-style update, and it's
  structurally weaker: naively *iterating* it (`u_n <- u_`, resolve, repeat)
  diverges exponentially by Re~2000 (`max|u|` observed going 20 → 86 → 3,683
  → ... → 1e117 over 8 iterations), whereas real Picard converges cleanly
  (observed: 6 iterations to `rel_change<1e-6` at Re=10). It also silently
  broke downstream consistency: a SUPG stabilization term built the same
  way (residual evaluated entirely at `u_n`, added only to the RHS) had no
  coupling to the actual unknown, so it could only ever perturb the
  solution — not stabilize it against the operator being solved — which is
  exactly why it visibly biased the velocity profile away from the
  (better-behaved) no-SUPG variant instead of just damping oscillations.
  General lesson: a "convection = frozen-at-`u_n`, forcing-term-only"
  scheme is a cheap trick worth trying to see how a solver breaks, but
  it is not a substitute for even the simplest consistent nonlinear
  iteration — and any *consistent* stabilization term added on top (SUPG,
  etc.) needs to be built from the same linearization the base equation
  uses, or it stops being consistent too.

- **[cavity]** — **Fixing a stabilization term's Picard-consistency (see the entry
  above) is necessary but not sufficient — a *consistent* stabilization term must
  also be built from the FULL equation residual, or it still biases the solution.**
  After making `fem_weak_SUPG.py`'s SUPG term properly bilinear-in-`u` (consistent
  with the fixed Picard scheme), it *still* undershoots both the Ghia et al. (1982)
  benchmark and the plain no-SUPG variant at every extremum, confirmed at Re=400 on
  the coarse mesh (`results/Re_400/coarse/dat/centerline.png`) — e.g. the `u`-vertical
  trough sits at `~-1.3` for Ghia/strong-BC, `~-1.25` for the no-SUPG weak variant, but
  only `~-0.5` for the SUPG one. This is the *opposite* of "SUPG's effect is too small
  to matter" (that was true at Re=100, and is expected/correct there — Pe_h<1
  everywhere, nothing to stabilize); at Re=400 SUPG is clearly engaged (46.6% of cells
  have local Peclet number `Pe_h>1`) and the effect is large and wrong-signed.
  **Root cause: `R1 = dot(u_n, nabla_grad(u))` only contains the convective part of the
  momentum residual** — the diffusive term (`-nu*div(grad(u))`) was deliberately
  omitted (differentiating a P2 field's second derivatives inside a bilinear form
  scaled by an unbounded `tau` severely ill-conditions the system), and the
  pressure-gradient term was never included either. The textbook guarantee that
  SUPG only damps genuine error — because a *true* residual is exactly zero at the
  exact solution, so the stabilization term vanishes there too — requires the FULL
  residual. A partial residual is not zero at the exact solution, so the "stabilization"
  term doesn't vanish there either; it adds a persistent, systematic offset instead of
  a correction that disappears as the solution converges. General lesson: before
  trusting a Petrov-Galerkin/residual-based stabilization scheme, check that every
  term of the actual PDE operator appears in the residual it's built from — leaving
  one out for numerical-conditioning reasons is a legitimate engineering tradeoff, but
  it silently converts "stabilizer" into "biased forcing term," and the two look
  identical from the diagnostics alone (smaller oscillations) without checking against
  an independent benchmark.
  **Fixed 2026-10-21, after comparing against `cylinder/fem/solver/solvers.py`'s
  `ChorinWeak` class**, which turns out to already do this correctly: its own SUPG
  residual *does* include the viscous term (`R1 = rho/k*(u-u_n) +
  rho*dot(u_n,nabla_grad(0.5*(u+u_n))) - 0.5*mu*div(grad(u+u_n)) - f` — unsteady +
  convection + viscous + forcing, i.e. the actual full residual of the equation that
  solver's Step 1 solves), and its `tau` has an extra `(2/dt)**2` term cavity's lacked
  entirely (cavity is steady, so there's no real `dt` anywhere). That term is what let
  cylinder get away with including the diffusive term at all: `(2/dt)**2` puts a hard
  ceiling on `tau` regardless of how small `h` or `nu` get, which is exactly the
  "unbounded tau" cavity's own comment flagged as the problem. Applied the same fix to
  `fem_weak_SUPG.py`: added the viscous term back into `R1`
  (`R1 = dot(u_n, nabla_grad(u)) - nu*div(grad(u))`, sign matching the strong-form
  momentum equation `-nu*div(grad(u)) + (u.grad)u + grad(p) = 0`), and introduced an
  artificial `dt_pseudo = L/U_lid` (domain-crossing time at the lid speed) purely to
  add the same `(2/dt_pseudo)**2` cap to `tau` — a number with no physical role
  anywhere else in this steady solver, existing solely to bound the stabilization
  parameter the way a real timestep incidentally does in an unsteady scheme. This is a
  general pattern worth remembering: **pseudo-transient continuation** (inventing an
  artificial timescale purely to regularize/bound a steady-state computation) is a
  legitimate, standard technique, not a hack — borrow the physically-motivated bound a
  transient scheme gets "for free" even when solving a genuinely steady problem. The
  pressure-gradient term `grad(p)` was still not included, and re-measuring (Re=400,
  coarse mesh) showed the undershoot barely changed: `strong`=`-1.246`, plain
  `weak`=`-0.969` (22% low), `weak_SUPG`=`-0.720` (42% low) — still clearly worse than
  the no-SUPG variant, by roughly the same margin as before this fix. So the
  diffusive-term incompleteness was *not* the dominant cause of the bias after all.
  Good reminder that fixing one identified incompleteness in a residual-based scheme
  doesn't let you assume the bias is gone — verify against the benchmark again,
  because another missing term can produce a very similar symptom.
  **Added `grad(p)` to `R1` too, same day, and this time it fixed it.** To be precise
  about why `grad(p)` belongs there: this has nothing to do with PSPG
  (pressure-stabilizing Petrov-Galerkin, which exists to fix the LBB/inf-sup
  instability of *equal-order* velocity-pressure spaces like P1-P1) — this solver uses
  Taylor-Hood (P2 velocity, P1 pressure), already LBB-stable, no PSPG needed. It's
  simply that SUPG's own `R1` is supposed to be the full momentum-equation residual,
  and `grad(p)` is part of that equation regardless of inf-sup considerations.
  Cylinder's `R1` skips it too, but legitimately — Chorin's splitting means Step 1's
  own equation genuinely has no pressure term in it. Cavity's monolithic
  velocity-pressure coupling has no such excuse. Mechanically, since `R1` now depends
  on *both* trial functions, the SUPG term no longer lives entirely in the
  velocity-velocity block: the `u`-dependent piece stays in `nitsche_a00`, and the new
  `grad(p)` piece becomes its own bilinear form in `(w, p)`, added to the existing
  `pressure` term in the velocity-pressure block instead
  (`a = [[diffusive + nitsche_a00 + convective, pressure + supg_pressure], [...]]`).
  Re-measured again (Re=400, coarse mesh): `weak_SUPG` moved from `-0.720` (42% low)
  to `-0.976` (22% low) — now matching plain `weak` (`-0.969`) almost exactly, on both
  the `u`-vertical and `v`-horizontal panels. The systematic undershoot is gone; the
  remaining `weak`/`weak_SUPG` vs. `strong`/Ghia gap is just the ordinary
  Nitsche-penalty approximation bias, not a SUPG artifact. Confirms the diagnosis
  cleanly: once a stabilization term is built from a genuinely complete residual, it
  behaves the way the textbook guarantee promises — agreeing closely with the
  no-stabilization result wherever the solution doesn't actually need stabilizing.

- **[cylinder]** — **MPI parallelism for FEniCSx/dolfinx scripts is set at the launcher,
  not in the script.** dolfinx scripts (e.g. `cylinder/fem/dokken/navStokes.py`) import
  `mpi4py.MPI` and pass `MPI.COMM_WORLD` into mesh reads/solves — that communicator adapts
  to however many ranks it's given, so there's no `n_cores` variable to edit anywhere in
  the script. Core count is chosen at invocation instead: `mpirun -np N python3 script.py`
  (or `mpiexec -n N python3 script.py`) — same pattern as OpenFOAM's `mpirun -np 16
  foamRun -parallel`. dolfinx then partitions the mesh across the N ranks automatically;
  no code change needed to go from 1 core to N. Two conventions that go with this, both
  already used in the dokken script: (1) guard any print/plot/file-write with `if
  mesh.comm.rank == 0:` — otherwise every rank does it redundantly, or N ranks race to
  write the same file; (2) to combine a per-rank scalar (e.g. a drag/lift integral
  computed over just that rank's mesh partition) into one value, use
  `mesh.comm.gather(value, root=0)` and sum/process the resulting list on rank 0 only.
