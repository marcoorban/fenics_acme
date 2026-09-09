# navStokes — Project Log

_Started: 2026-09-09_

This is the **project-wide** log: tooling, workflow, and technical
gotchas that apply across cases, or that were learned in one case but
generalize to any future one. Case-specific decisions/tasks (mesh
parameters, boundary conditions, Reynolds numbers, benchmark comparisons
for that specific geometry) stay in the per-case logs:

- `cavity/log.md` — lid-driven cavity case.
- `cylinder/log.md` — flow past a cylinder (not yet started).
- `cavity/plot_notes.md` — detailed rationale for the `plots.gp` styling
  pass (kept separate since it's long and cavity-specific in its
  examples, even though the *convention* it establishes is project-wide
  — synthesized into `#plotting` below).

## How entries here are tagged

Each entry has a **scope** tag and one or more **topic** tags:

- **Scope** — `[project]` (repo/tooling-wide, not about any one case) or
  `[cavity]`/`[cylinder]` (the lesson was learned working on that case,
  but is written up here because it's expected to generalize to other
  cases too).
- **Topic** — `#openFoam`, `#gmsh`, `#git`, `#docker`, `#plotting`,
  `#fem` — the subsystem the note is actually about.

## #docker

- **[project]** — `oforg.fish` (the Docker wrapper for
  `microfluidica/openfoam:org`) must run as `--user (id -u):(id -g) -e
  HOME=/tmp`, not root. `-v $PWD:/project` bind mounts share the host's
  user namespace, so a container run without `--user` writes
  root-owned output back to the host; `HOME=/tmp` is needed alongside it
  because the host UID has no entry in the container's `/etc/passwd`.
- **[project]** — "permission denied" talking to
  `/var/run/docker.sock` after adding yourself to the `docker` group
  usually means the *current shell session* hasn't picked up the new
  group membership yet (group lists are fixed at login/process
  creation) — confirm with `id`/`groups` showing `docker` missing, then
  fix with `newgrp docker` (immediate, current shell only) or a fresh
  login/terminal (persists). Run in your own terminal, not through a
  tool that spawns isolated one-off shells per command — the group
  change won't carry over between separate invocations there.

## #git

- **[project]** — `.gitignore`'s OpenFOAM section is anchored on the
  path *segment name* `openFoam` at any depth
  (`navStokes/**/openFoam/**/...`), not a literal prefix
  (`navStokes/openFoam/**/...`). A literal-prefix anchor breaks silently
  (no error, previously-ignored files just become trackable again) the
  moment a case's directory layout changes — which happened once already
  (case moved from `navStokes/openFoam/<case>/` to
  `navStokes/<case>/openFoam/`, ~653 generated files went from ignored
  to staged before the fix).
- **[project]** — `git status`'s `R` (rename) lines are a *content
  similarity* heuristic, not proof of an actual `mv`/`git mv`. Seen
  pairing an unrelated staged deletion with an unrelated staged addition
  purely because both files were short and similar. Cross-check with
  `find`/`git log` on the real paths before trusting an `R` line.
- **[project]** — Untracked files don't travel with a git-aware
  directory reorg — they're invisible to `git mv`/rename-pairing, so
  they silently get left behind at the old path while every tracked
  sibling moves cleanly. Check for untracked files at the old location
  specifically after any large restructuring.
- **[project]** — A delete/modify merge conflict (`deleted by us`) means
  your branch removed a file while the other branch kept/changed it —
  git checks "their" version into the working tree with no conflict
  markers (nothing to textually merge), but still requires an explicit
  `git add`/`git rm` to resolve. If the deletion was intentional,
  `git rm <path>` confirms it and stages the resolution.

## #openFoam

- **[project]** — `fvSolution`'s `solvers{}` `tolerance`/`relTol` and
  `PIMPLE{}`/`SIMPLE{}`'s `residualControl` are **not** the same kind of
  "convergence" setting, easy as that is to assume:
  - `tolerance`/`relTol` (per solver, e.g. `p`, `U`) control how tightly
    that *one linear system* is solved **within a single time step** —
    numerical accuracy of that step's solve, nothing to do with whether
    the physical solution has stopped evolving.
  - `residualControl` (in `PIMPLE`/`SIMPLE`) is the actual steady-state
    detector: it checks each listed field's *initial* residual every
    time step and auto-stops the run once all of them are below
    threshold simultaneously — the right tool for "stop once it's
    converged" instead of always running to `endTime`.
- **[cavity]** — `gmshToFoam` infers a boundary's OpenFOAM `type`
  (`wall` vs. generic `patch`) heuristically from the *physical-surface
  name* — a group literally named `"wall"` gets `type wall`; one named
  `"lid"` (despite being a real, moving wall) gets `type patch`.
  Harmless for a laminar Dirichlet velocity BC, but would matter for
  wall-function/`wallShearStress`/`yPlus` post-processing on that patch.
  Worth checking on any future case's `constant/polyMesh/boundary`.
- **[cavity]** — `internalField` in a field file (`U`, `p`, ...) is
  per-*cell*, in the mesh's own internal numbering — not a
  reshape-able `i,j` grid, even on a structured hex mesh. For a profile
  along `x=const`/`y=const`, use the `sample` function object
  (`system/sampleDict`, `type sets;`, run via `postProcess -func
  sampleDict`) — it interpolates via the mesh's actual connectivity
  (`interpolationScheme cellPoint`) rather than assuming a scan order.

## #gmsh

- **[cavity]** — `gmsh.model.geo.extrude(...)` (Python API) returns
  entities as `[far cap, volume, lateral surfaces in Curve Loop
  order]` — same convention as the `.geo`-script `Extrude{}` command.
  Verify per-geometry against the resulting `constant/polyMesh/boundary`
  face counts rather than trusting the mapping from memory — renumbering
  points/lines changes it.

## #plotting

Full rationale and options considered live in `cavity/plot_notes.md`;
synthesized here as the **project's plotting convention** (applies to
any future case's validation plots, e.g. cylinder drag/lift/Strouhal
comparisons):

- **[project]** — Validation-figure convention: tabulated *benchmark*
  data (few discrete points) as markers only; the *current simulation*
  (densely sampled) as a line only. Both `with line` misleadingly
  implies a coarse reference table is a continuous measurement.
- **[project]** — Okabe-Ito colorblind-safe palette (black markers +
  `#0072B2` blue line) instead of gnuplot's default red/green cycling.
- **[project]** — `pngcairo enhanced font "Helvetica,14"` +
  italicized axis labels (`{/Italic y}`) — matches ISO 80000 /
  journal-style convention for physical-quantity symbols. Vector output
  (`pdfcairo`/`epslatex`) considered for a genuinely journal-ready figure
  but not applied yet — raster PNG kept for now.
- **[cavity]** — Before plotting a *raw* OpenFOAM field value against a
  *published benchmark* table, check whether the benchmark is
  non-dimensionalized. Ghia et al. (1982)'s `u`/`v` are scaled by the
  lid speed (`u* = u/U_lid`) — invisible when `U_lid=1` (numerically
  identical to the raw value), but a stark mismatch once `U_lid` changes
  (e.g. Re=400 run uses `U_lid=4`: raw `U_x` ranges to 4.0, Ghia's
  column ranges to 1.0 — same shape, 4× the scale). Verified the
  normalization itself directly from the data: the lid-boundary row
  (`y=1`) is exactly `1.00000` across *every* Reynolds-number column in
  Ghia's table, which only holds if velocities are reported as a
  fraction of the (fixed, non-dimensional) lid speed. General lesson:
  check a benchmark table's boundary-condition row for this kind of
  tell before assuming raw and tabulated values are directly comparable.
