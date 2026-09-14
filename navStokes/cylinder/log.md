# navStokes/cylinder — Log

_Started: 2026-09-10_

## Action Items

- [x] **✅ 2026-09-14 — DFG 2D-3 (Re=100, time-ramped inlet) run, analyzed, matches the published benchmark
      closely.** Same case/mesh as `2D-2`, `0/U`'s inlet BC extended with a `sin(pi*t/8)` time factor
      (`Um=1.5` unchanged), `controlDict` `endTime` `5 -> 8s`. `Cd max=2.9632` at `t=3.92s` (ref. `~2.9483`,
      ~0.5% off), `Cl max=0.4803` at `t=5.72s` (ref. `~0.4651`, ~3% off), `delta_p(t=8)=-0.131359`. Full
      numbers/method in `results.csv`. Snapshot (mesh + `t=3.92s`/`5.72s`/`8s` field dirs + postprocessing
      data) in `openFoam-results/2D-3_unsteady/`; `plots/forceCoeffs.png` (via `plot2d-3.gp`, `Cd`/`Cl` vs.
      `t` with horizontal `9a`-reference dashed lines) also copied there and is the only thing from that
      directory actually git-tracked (see the `openFoam-results/` gitignore correction below).

- [x] **✅ 2026-09-14 — DFG 2D-2 (Re=100, unsteady) analyzed and matches the published benchmark closely**
      (`St=0.3007` vs. ref. `~0.300`, `Cd mean=3.202` vs. ref. `~3.2`, `Cl max=1.001` vs. ref. `~1.0`).
      `plot_results.py` already accounts for the amplitude-growth transient (see Learnings) by using a
      `t>=4.0s` saturated window for the extrema/`delta_p` stats and a wider `t>=3.0s` window for Strouhal
      (frequency is unaffected by the amplitude still growing). Not a loose end — if a tighter number is
      ever wanted, extend `endTime` past `5s`, but current values are already a close match.

- [x] **✅ 2026-09-14 — RESOLVED, via mesh regrading rather than the corrector fix in isolation.** Rather
      than retesting `nNonOrthogonalCorrectors` alone at `deltaT~1e-5` as planned below, the mesh itself
      was regraded (`Distance`+`Threshold` field, `recombine=True`) — see Learnings. **This means the old
      question ("does the corrector fix alone work on the original fine mesh?") was never actually answered
      and no longer will be** — the fine, ungraded mesh this bug was originally diagnosed on has been
      superseded. Not a loose end worth chasing further: the regraded mesh is better on every axis (quality,
      cell count, now confirmed stability), so there's no reason to go back and isolate the old variable.

- [x] `deltaT` is back at `1e-5`, confirmed stable and converging (`residualControl` triggers around
      `t≈0.3-0.4s` on the coarse meshes tested so far).

- [ ] `system/changeDictionaryDict` + the `changeDictionary` utility (used 2026-09-14 to force `empty`/
      `wall` patch types after `gmshToFoam`) is now **REMOVED — do not use `changeDictionary` again, it's
      deprecated.** Superseded by `openFoam/mesh.sh` (user-authored): runs `gmshToFoam` then fixes the same
      patch types via `foamDictionary -set`. Always run mesh conversion through `mesh.sh` from now on, not
      manual `gmshToFoam` calls.

## Objective

OpenFOAM case for classic 2D flow past a cylinder — the DFG 2D-1 (Schäfer–Turek) benchmark geometry: channel
`2.2 x 0.41`, cylinder diameter `D=0.1` at `(0.2, 0.2)`, asymmetric gaps to the walls (`0.15` below, `0.16`
above) to break vertical symmetry and speed up Kármán vortex shedding. Mesh built with `mesh/mesh.py` (gmsh
Python API, unstructured triangles refined near the cylinder, thin single-layer extrusion for OpenFOAM's
pseudo-2D convention, same pattern as `cavity/mesh/mesh.py`), converted via `gmshToFoam`. Run with `foamRun
-solver incompressibleFluid`, `simulationType laminar`.

## Key decisions

- **2026-09-10** — `constant/physicalProperties`'s `nu` was `1e-2` (stale, copied from `cavity`). Updated to
  `1e-3`. Why: for a parabolic profile the Reynolds number convention uses the *mean* velocity, not the
  peak — `Ubar = (2/3) Um`. With `Um=0.3`, `D=0.1`: `Ubar=0.2`, so `Re = Ubar*D/nu = 20` requires `nu =
  1e-3`, not `1e-2`.

- **2026-09-10** — `0/U` and `0/p` were stale copies from the cavity case (`lid` patch, which doesn't exist
  on this mesh — `inlet`, `outlet`, `cylinder` were entirely undefined). Rewrote both against this mesh's
  actual patches:
  - `inlet`: `codedFixedValue` applying the parabolic (Poiseuille) profile pointwise, `U(0,y) = 4 Um y (H-y)
    / H^2`, `V=0`. `Um` and `H` are named constants at the top of the `code` block — `Um` is the single
    value to edit for different flow conditions; OpenFOAM hashes the code string and auto-recompiles the BC
    when it changes.
  - `outlet`: `U` `zeroGradient`, `p` `fixedValue 0` (open-boundary reference pressure).
  - `wall`/`cylinder`: `noSlip` for `U`, `zeroGradient` for `p`.
  - `empty`: unchanged (front/back caps).

- **Mesh**: unstructured (not structured/transfinite like `cavity`) — the circular cutout makes a clean
  structured grid impractical. Physical groups: `inlet`/`outlet`/`wall`(top+bottom)/`cylinder`/`empty`
  (front-back caps)/`domain`, matching `cavity`'s `lid`/`wall`/`empty`/`domain` naming convention.

## Learnings (technical gotchas worth not re-deriving)

- **🚨 2026-09-14 — `reconstructPar -latestTime` only reconstructs the single latest time, and
  `foamPostProcess ... -latestTime` (run in serial, against the reconstructed case) can only see whatever
  times actually got reconstructed** — for the `2D-3` run this meant `forceCoeffsDict` initially produced a
  single-row `.dat` file (just `t=8`), not the full time series needed to find `Cd`/`Cl` maxima over the
  run. The `-latestTime` flag on `foamPostProcess` does NOT mean "sweep every written time and use the
  latest mesh/field state for setup" the way it might read — it means "only evaluate at whichever time is
  latest among what's actually present on disk". **Fix, and the better approach in general**: the full
  per-time-step field data was never lost — it's still in the `processorN/` parallel-decomposed directories
  (`foamRun`'s native output before any reconstruction). Ran `foamPostProcess -func forceCoeffsDict
  -parallel` (via `mpirun`) directly against that decomposed data instead — computes the same force
  integrals over each subdomain's slice of the `cylinder` patch and combines via MPI reduction, so it gets
  the complete `402`-point time series without ever reconstructing full field snapshots for every step
  (which would've been `~1.9G` at `4.7M`/step here, entirely wasted just to compute two scalars per step).
  **General lesson: for scalar/integral postprocessing (`forceCoeffs`, `pressureDrop`, etc.), prefer
  `-parallel` against the decomposed case over reconstructing first** — reconstruction is only actually
  needed when you want full field snapshots (for visualization, or to save to `openFoam-results/`), not for
  postprocessing utilities that reduce a field to a handful of numbers per time. Needed a *targeted*
  `reconstructPar -time 3.92,5.72` (comma list, not `-latestTime`) afterward anyway, to get the actual
  `U`/`p` field snapshots at the `Cd`/`Cl`-max times for `openFoam-results/`.

- **2026-09-14** — `openFoam-results/` is where a finished case's *permanent, git-tracked* artifacts go,
  separate from `openFoam/` (whose raw solver output is entirely gitignored and regenerates on every
  rerun/mesh regen). First pass saved raw solver output here too (mesh + full time-directory field data) —
  for `2D-2_unsteady` (Re=100, periodic) that meant the mesh + **every written time directory spanning two
  full shedding periods** (`t=4.32` to `t=5s`, 34 directories at the existing `0.02s` write spacing, bounded
  by the same `Cl`-maximum times used in the `delta_p` calc above), enough to reconstruct/animate a full
  cycle later — but this ballooned to `176M` (mesh `13M` + `35` field snapshots at `4.7M` each), an order of
  magnitude more than `2D-1_steady`'s single-snapshot `14M`, and isn't something worth tracking in git.
  **Corrected**: `openFoam-results/*/[0-9]*/` (time directories) and `openFoam-results/*/postProcessing/`
  added to `.gitignore` (mesh/dat files were already covered by existing blanket `*.msh`/`*.dat` rules) —
  careful to anchor on `*/[0-9]*/` (exactly one directory down from `openFoam-results/`), not `**/[0-9]*/`,
  since the case directory names themselves (`2D-2_unsteady`) also start with a digit and would otherwise
  match too, wiping out the whole case directory instead of just its time subdirectories. Raw data stays on
  local disk (not deleted, just untracked) for now, in case it's useful for a local animation later. What
  actually stays git-tracked in `openFoam-results/` going forward is small summary artifacts only — for
  `2D-2_unsteady`, the three validation plots (`cd_vs_t.png`/`cl_vs_t.png`/`dp_vs_t.png`, copied in from
  `cylinder/plots/`) rather than the data that produced them. `2D-1_steady`'s single time-directory snapshot
  is now also gitignored by the same rule — its useful content (the converged `Cd`/`Cl`/`delta_p` numbers)
  already lives in `results.csv`, not in the raw field data, so nothing of value was lost by no longer
  tracking it. Worth keeping "small derived summary tracked, bulky raw solver output gitignored" as the
  pattern for `2D-3` (oscillating inlet, unsteady, not started) too, rather than repeating the raw-data-
  first approach here.

- **2026-09-14** — Full time history (not just the periodic tail) shows the classic unsteady-cylinder
  startup transient clearly: `Cd` spikes to `~2.5` in the first couple of timesteps then relaxes to `~3.0`
  by `t~0.2s` and creeps slowly up toward its limit-cycle mean from there; `Cl` stays near-zero (symmetric
  wake) until the wake instability grows enough to break symmetry, then the shedding oscillation grows
  quasi-exponentially in amplitude until saturating into the limit cycle around `t~3-4s`. Worth keeping in
  mind for any future case: window any periodicity/FFT-based analysis (Strouhal number, amplitude) well
  after this initial transient, not from `t=0`.

- **2026-09-14** — `forceCoeffsDict`/`pressureDropDict` postprocessing (`foamPostProcess -func ...
  -latestTime` in `foamRun.fish`, run after `reconstructPar`) takes noticeably longer than the solve step's
  final `ExecutionTime` would suggest — it re-walks *every* written time directory (251 of them here,
  `writeInterval=200` at `deltaT=1e-4`) to build the full `.dat` time series, not just the one named by
  `-latestTime` (that flag only affects which *mesh/field state* is used to construct the function object,
  not which times it evaluates over). Don't assume a `.dat` file is complete just because the solve step's
  `log.foamRun` already shows `End`/`Finalising parallel run` — check the `.dat` file's own last time value
  against the expected `endTime`, or that its row count matches the number of time directories on disk,
  before treating post-processed results as final.

- **🚨 2026-09-14 — `delta_p` was initially computed wrong: a plain time-average over the saturated window,
  not the actual DFG benchmark definition.** The benchmark defines `delta_p = p(front) - p(rear)` at a
  *specific instant*, `t0 + T/2`, where `t0` is a time at which `Cl` attains a local maximum — not a
  statistic over an arbitrary window. A plain mean/range over `t>=4.0s` mixes pressure-drop values from all
  phases of the shedding cycle together (mean `2.4530`, range `[2.4115, 2.4946]` — the range alone shows how
  much this varies over a cycle, so a window-average was never going to be a meaningful number to compare
  against the benchmark). **Fixed** in `plot_results.py` (`delta_p_at_half_period`): find each `Cl` local
  max in the saturated region (parabolic-vertex refinement across 3 samples for sub-sample-accurate peak
  timing, since `deltaT=1e-4` sampling alone isn't fine enough to land exactly on the true peak), evaluate
  `dp` at `t0 + T/2` via linear interpolation into the pressure-probe series (`T` from the Strouhal calc
  above), and average across every `t0` whose `t0+T/2` still falls inside the recorded series (2 of the 3
  `Cl` maxima found at `t0>=4.0s` qualify; the third's half-period point lands past `t=5s`). Result:
  `delta_p=2.469033`, samples `[2.468752, 2.469314]` — spread `<0.001`, i.e. essentially reproducible
  cycle-to-cycle, unlike the window-average approach. `results.csv` corrected accordingly. **General lesson:
  a benchmark quantity defined at a specific phase of a periodic signal is not the same as that quantity's
  average over the cycle** — check the actual definition before reaching for a window mean/min/max,
  especially for anything that (unlike `Cd`/`Cl` extrema, which genuinely *are* per-cycle extrema) is
  defined relative to another signal's phase rather than as a standalone statistic.

- **2026-09-14** — `plot_results.py` uses **two different saturation windows on purpose**: `T_SATURATED=
  4.0s` (tight) for `Cd`/`Cl` extrema and `delta_p`, since those are sensitive to the shedding amplitude
  still growing from rest early on; `T_STROUHAL=3.0s` (wider) for the period-averaging, since the shedding
  *frequency* doesn't change while the amplitude is still growing — more zero-crossings in the wider window
  means a better-averaged period estimate, with no accuracy cost from including the "not fully saturated
  yet" portion. Worth keeping this asymmetric-window pattern in mind for any future unsteady/periodic case
  analysis: don't apply one blanket "steady-state cutoff" time to every derived quantity if some are
  amplitude-sensitive and others aren't.

- **2026-09-14 — Re=100 result: excellent match to the published DFG 2D-2 reference ranges.**
  `plot_results.py` computes: `Strouhal St=0.3007` (ref. `~0.300`) via positive-going `Cl` zero-crossings
  over `t>=3.0s` (5 periods, mean `T=0.3325s`); over the saturated `t>=4.0s` tail, `Cd` mean `3.2022`, range
  `[3.1711, 3.2337]` (ref. `c_D,mean~3.2`), `Cl` range `[-1.0014, 0.9675]` (ref. `c_L,max~1.0`). `delta_p`
  is handled separately below (2026-09-14 correction entry) — it's not a window statistic like the others.
  Full numbers/method in `results.csv`. Two sets of plots: `plots/re100_force_coefficients.png`
  (`plot_results.py`, full `0-5s` run, shows the startup transient) and
  `plots/cd_vs_t.png`/`cl_vs_t.png`/`dp_vs_t.png` (`plots.gp`, saturated-window zoom, project `#plotting`
  convention — horizontal dashed reference lines instead of benchmark markers here, since DFG 2D-2 publishes
  only scalar `St`/`c_D`/`c_L` reference *ranges*, not a tabulated time series to overlay the way `cavity`'s
  Ghia profiles are).

- **2026-09-14 — moved on to the DFG 2D-2 benchmark (Re=100, unsteady vortex shedding)**, after `2D-1`
  (Re=20, steady) converged cleanly. Only `Um` (`0/U`, `0.3 -> 1.5`) and `nu`
  (`constant/physicalProperties`, matched to keep `D=0.1` fixed) needed to change to retarget `Re`; also
  bumped `forceCoeffsDict`'s `magUInf` to the new `Ubar=1.0` (see the reference-value comment already in
  that file — easy to forget this when only `Um` is changed, since `Ubar=(2/3)Um` is a derived quantity, not
  read from `0/U` automatically). Mesh regenerated finer than the `2D-1` run (`size_far=0.0097`,
  `size_cyl=0.00194` -> 80195 cells) to resolve the shed vortices adequately; `endTime` set to a fixed `5s`
  since `PIMPLE.residualControl` is meaningless here by design — the flow never reaches a steady state, so
  it (correctly) never triggers and the run always goes to `endTime`.

- First two mesh-sizing iterations toward a 50k-75k cell target overshot low: naive `1/size²` scaling from a
  known `(size, cells)` data point undershot badly (predicted ~59k cells at `size_far=0.02, size_cyl=
  0.004`, actual was 19075) — the far-field region's contribution didn't scale as cleanly as assumed,
  likely due to the graded-transition zone (`dist_min`→`dist_max`) taking up a fixed-size chunk of the
  domain regardless of `size_far`/`size_cyl`. Needed a second empirical correction (scale down by
  `sqrt(target/actual)` from the actual measured result, not the original estimate) to land in range
  (`size_far=0.011, size_cyl=0.0022` → 62373 cells). **Treat mesh-size-to-cell-count scaling as empirical,
  not analytically predictable, once a background field with a fixed transition distance is involved** —
  iterate from the actual `checkMesh` result, not from a formula.

- **gmsh's own `"N nodes M elements"` summary line is not a proxy for the actual OpenFOAM cell count** — it
  bundles 2D surface/cap triangles and boundary-patch faces together with the 3D cells. Seen ranging from
  `cells ≈ 0.33×` to `≈ 0.95×` of that total depending on resolution — no fixed ratio to convert by. **Always
  check `checkMesh`'s own `cells:` line**, not gmsh's summary, when targeting a specific cell count.

- **`foamPostProcess -func forceCoeffsDict` failed with `Could not find U, p` even though those files
  genuinely existed on disk** at the requested time directory (verified: valid FoamFile headers, correct
  `volScalarField`/`volVectorField` class, real data). Fix: add `-solver incompressibleFluid` to the
  `foamPostProcess` call — apparently this utility needs that hint to construct the solver's field registry
  before a function object can look fields up in it; without it, the mesh loads fine but the fields are
  never registered. Updated in `foamRun.fish`.

- **`gmshToFoam`'s automatic `wall`/`empty` patch-type detection is NOT reliable — don't trust it, verify
  `constant/polyMesh/boundary` after every conversion.** It worked on the original fine mesh (and still
  works on `cavity`), but silently fell back to generic `type patch;` for *every* patch — including
  `empty` — on the first regraded/coarser mesh, with no error or warning at conversion time. This is fatal
  downstream, not benign: an `empty`-type field BC requires the *mesh* patch to be constraint-type `empty`,
  so `decomposePar` failed with `FOAM FATAL IO ERROR: patch type 'patch' not constraint type 'empty'`. Root
  cause not fully pinned down (tried `recombine=True` as a fix on the hypothesis that mismatched front/back
  triangulation from tet-subdivision was confusing the detection — quality improved but the type bug
  persisted, so that wasn't it either). **Fixed by forcing the types explicitly** rather than depending on
  the heuristic at all: first via `system/changeDictionaryDict` + `changeDictionary` (now REMOVED —
  deprecated, see Action Items), now via the user's `openFoam/mesh.sh` (`gmshToFoam` then `foamDictionary
  -set` on `constant/polyMesh/boundary`). **Always verify patch types after any mesh regeneration** — this
  isn't a one-time fix, it's a per-regeneration risk.

- **🚨 2026-09-10/14 — Courant number "blowing up" on the early meshes was never a timestep problem, and the
  real fix ended up being the mesh itself, not the `nNonOrthogonalCorrectors` tweak tried first.** Two
  distinct blow-ups showed up, told apart by checking **mean vs. max Courant number separately** (the
  general lesson worth keeping): a domain-wide-high *mean* is a genuine CFL/timestep issue; a *mean* pinned
  near-zero with an exploding *max* is a localized numerical-instability signature that no amount of
  shrinking `deltaT` fixes — confirmed directly, dropping `deltaT` `1e-5→1e-6` still blew up at the same
  physical time. The first (mean-near-zero) blow-up traced to `checkMesh` non-orthogonality (`max 83.9°`,
  ~12% severely non-orthogonal faces) with `nNonOrthogonalCorrectors=0`; bumping that to `2` was tried, but
  a follow-up test bundled it with an unrelated 100× `deltaT` jump in the same run, so it was never cleanly
  confirmed in isolation (that second blow-up had an elevated *mean* too — a plain CFL issue, a different
  failure mode entirely). **What actually fixed it: regrading the mesh**, not the corrector count. Replaced
  the old sharp two-value point sizing (jumping `0.01→0.001` right at the cylinder surface) with a
  `Distance` field on the cylinder arcs feeding a `Threshold` field — smooth linear grading from `size_cyl`
  at the surface out to `size_far` beyond `dist_max=3D` — plus `recombine=True` on the extrude (missing
  from `cylinder/mesh.py`, unlike `cavity`'s, so cylinder had been producing tet-subdivided cells instead of
  clean prisms). Together these dropped non-orthogonality to max `21.9°`/avg `4.2°` with zero severely
  non-orthogonal faces. Validated the new grading on a deliberately coarse mesh first (confirmed no
  blow-up), then refined the cell count for the real run.

- **2026-09-10** — `Allclean` never cleaned `processorN/` directories, only top-level time dirs!! After
  several `mesh/mesh.py` revisions (domain rescaled to the DFG `2.2 x 0.41` geometry, `D=0.1`, mesh refined
  near the cylinder), `processor0-7` still held a full old solved run from an earlier, much coarser mesh
  (~2080 cells/subdomain vs. the current mesh's ~13352/subdomain across 8 subdomains). `reconstructPar
  -latestTime` grabbed that stale leftover time directory (the fresh run hadn't gotten anywhere near its
  `t≈2.77`) and tried to map its old-sized field data onto the newly-rebuilt decomposition addressing:
  `FOAM FATAL IO ERROR: size 2080 is not equal to the given value of 13352`!! Fixed by adding `rm -rf
  processor*/` to `Allclean` (see `[cylinder]`/`#openFoam` in the project log — same latent bug found in
  `cavity/Allclean` too and fixed there as well).

- Copying `0/` field files from another case (`cavity` → `cylinder`) carries over that case's patch names
  silently — no error until something tries to use the missing/extra patches (`decomposePar`, or just
  running with an incomplete `boundaryField`). Worth diffing new-case `0/` files against
  `constant/polyMesh/boundary`'s actual patch list before running, not just against the old case's files.

- For a parabolic/Poiseuille inlet profile, `Re` is conventionally reported against the *mean* velocity
  `Ubar = (2/3) U_peak`, not the centerline peak value used directly in the profile formula — easy to
  conflate since the profile formula itself is written in terms of the peak. Cross-check `nu` against the
  target `Re` using `Ubar`, not `Um`.

- `codedFixedValue` is the right tool for a spatially-varying boundary value (can't express `U(y)` as a
  single `fixedValue` vector) — code lives directly in the field file, compiled at runtime and auto-rebuilt
  when the code string's hash changes. Cheap way to keep one flow-condition parameter (`Um`) easily
  sweepable without a separate compiled library.

## References

- Schäfer, M., Turek, S. (1996), *Benchmark Computations of Laminar Flow Around a Cylinder*, in: Hirschel
  E.H. (eds) Flow Simulation with High-Performance Computers II. Notes on Numerical Fluid Mechanics, vol 52.
  Vieweg+Teubner Verlag. (DFG 2D-1/2-2/2-3 benchmark suite — this case targets 2D-1, `Re=20`, steady.)
