# Cavity solver scheme

All three solvers (`fem_strong.py`, `fem_weak_SUPG.py`, `fem_weak.py`) share the same
nonlinear-iteration scheme for the convective term. Documented once here instead of
repeated in every script's comments.

## Linearization scheme

**Picard iteration** (successive substitution). The convective term `(u.grad)u` is
nonlinear because the unknown `u` plays both roles at once -- the transported quantity
and the transport velocity. Picard breaks that by replacing one copy of `u` with a
known/lagged field `u_n` (so `u_n` plays the role of a prescribed transport velocity,
exactly like `a` in a passive advection-diffusion problem), keeping the other copy as
the real trial function `u`. Each individual linear solve is then genuinely linear in
`u`; solving, setting `u_n <- u`, and repeating drives `u_n` toward `u`, and at that
fixed point the Picard-linearized equation coincides with the true nonlinear steady
Navier-Stokes equation.

Full derivation, and why the earlier "both copies of u frozen at u_n, moved entirely to
the RHS" scheme was *not* actually Picard (and why that mattered for SUPG consistency),
is in `navStokes/log.md`'s `#fem` entry, 2026-10-21.

## Reynolds-number continuation

Plain Picard started from rest (`u_n=0`) stalls in a low-amplitude oscillation at
Re=400 instead of converging -- confirmed empirically: 200 iterations, `rel_change`
plateaued around 2e-5, never reached tolerance. Fix: reach the target Re through a
geometric ramp instead of one jump from rest. Starting at `min(50, Re)` and doubling
each stage up to the target, each stage's Picard loop starts from the *previous*
stage's converged `u_n` rather than from rest. Easy cases (Re<=50) reduce to a single
direct stage, so this adds no overhead there.

## Parameters

Set at the top of each script's solver block (identical across all three scripts):

| Parameter | Value | Meaning |
|---|---|---|
| `picard_tol` | `5e-6` | relative L2 change in `u` between iterations that counts as converged (`\|\|u_new - u_n\|\| / \|\|u_new\|\|`) |
| `picard_max_iter` | `200` | cap per continuation stage, not per whole run |
| `picard_omega` | `0.5` | under-relaxation factor on the `u_n` update (`u_n <- u_n + omega*(u_new - u_n)`); `1.0` would be unrelaxed Picard, which stalls at higher Re (see above) |

`picard_tol` was originally `1e-6`; relaxed to `5e-6` after confirming (Re=400, normal
mesh) that the continuation+relaxation scheme gets to `rel_change` in the `3e-6`-`6e-6`
range and stays there -- `5e-6` is reachable, `1e-6` was not within a reasonable
iteration budget.

## Logging

Every solver run appends one line to `results/dat/solver_timing.dat` -- this one file
stays at the top level (not per-mesh) since it's designed to compare across every
method/mesh/Re combination in one place:

```
# mesh method Re iterations time_s omega converged
```

`iterations` and `time_s` are totals across *all* continuation stages for that run, not
just the final target-Re stage. `converged` is `1` if the final stage's `rel_change` was
under `picard_tol` when the loop exited, `0` if it hit `picard_max_iter` first.

## Results layout

```
results/Re_<re>/<mesh_label>/{strong,weak,weak_SUPG}/   -- VTK output (u.pvd, p.pvd, ...)
results/Re_<re>/<mesh_label>/dat/                        -- centerline.dat/_vertical.dat/.png
results/dat/solver_timing.dat                            -- shared across every Re/mesh (see above)
```

Reynolds number is the outermost grouping, then mesh, then BC method -- a run is
identified by `(Re, mesh, method)` in that order of how often each one changes across
a sweep (least often to most).

`<mesh_label>` comes from the `CAVITY_MESH` environment variable -- set once by each
`run_sweep_*.fish` script (`set -gx CAVITY_MESH <label>`) and read by every solver
script, `centerline_profile.py`, and `plot_centerline.gp` (via `system('echo -n
$CAVITY_MESH')`, since gnuplot has no `getenv()` in this version). Set it yourself
before running a script by hand: `set -gx CAVITY_MESH coarse` (fish) or `export
CAVITY_MESH=coarse` (bash), and every script agrees on the same results folder without
editing anything. If `CAVITY_MESH` is unset, the Python scripts fall back to guessing
from `mesh.py`'s cell-count filename (see `MESH_LABELS` in each), and
`plot_centerline.gp` falls back to `normal`:

| cells | label | `mesh.py --scale` |
|---|---|---|
| 1968 | `coarser` | 2.514 |
| 2528 | `coarse` | 2.236 (`sqrt(5)`) |
| 5432 | `normal` | 1.414 (`sqrt(2)`) |
| 10504 | `fine` | 1.0 (baseline) |
| 15152 | `finest` | 0.816 (`1/sqrt(1.5)`) |

The three BC-method folders match the script names directly now: `fem_strong.py` ->
`strong`, `fem_weak_SUPG.py` (Nitsche + SUPG) -> `weak_SUPG`, `fem_weak.py` (Nitsche
only) -> `weak`.
