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

## SUPG residual completeness (fixed 2026-10-21)

`fem_weak_SUPG.py`'s SUPG term used to be built from an **incomplete residual** --
`R1 = dot(u_n, nabla_grad(u))` was the convective term only, with the diffusive term
(`-nu*div(grad(u))`) deliberately left out: including it meant differentiating the P2
velocity field's second derivatives inside a bilinear form scaled by an unbounded
`tau`, which severely ill-conditioned the system.

A *correctly consistent* stabilization term vanishes identically when evaluated at the
exact solution, since the true PDE residual is zero there by definition -- that's what
guarantees it only damps genuine error without biasing the answer. A partial residual
breaks that property, so it doesn't just damp oscillations, it injects a systematic
offset. **Confirmed empirically** (Re=400, coarse mesh): `weak_SUPG` undershot both the
Ghia et al. (1982) benchmark and the plain `weak` (no-SUPG) result at every extremum,
even in a regime where SUPG should matter a lot (46.6% of cells have local Peclet
number `Pe_h > 1` there, from `peclet_field.py`).

**Fix**: comparing against `cylinder/fem/solver/solvers.py`'s `ChorinWeak` class showed
its own SUPG residual already includes the viscous term, and its `tau` has an extra
`(2/dt)**2` term cavity's lacked -- that's what lets it get away with the diffusive
term at all, since `(2/dt)**2` caps `tau` regardless of how small `h`/`nu` get. Cavity
is steady (no real `dt`), so `dt_pseudo = L/U_lid` (domain-crossing time at the lid
speed) was introduced purely to provide the same cap -- a **pseudo-transient
continuation** device with no physical role anywhere else in this solver, existing
solely to bound `tau`. `R1` is now `dot(u_n, nabla_grad(u)) - nu*div(grad(u))`.

**Re-measured after that fix (Re=400, coarse mesh), and the undershoot was still
there**: at the `u`-vertical trough, `strong`=`-1.246`, plain `weak`=`-0.969` (22%
low), `weak_SUPG`=`-0.720` (42% low) -- `weak_SUPG` was still clearly worse than the
no-SUPG variant. So the diffusive-term incompleteness wasn't the dominant cause of the
bias; the remaining gap was `R1` still missing the pressure-gradient term `grad(p)`.
This has nothing to do with PSPG (pressure-stabilizing Petrov-Galerkin, which fixes
the LBB/inf-sup instability of *equal-order* velocity-pressure spaces like P1-P1) --
this solver uses Taylor-Hood (P2 velocity, P1 pressure), which already satisfies LBB
with no stabilization needed for that. The gap was narrower and simpler: SUPG's own
`R1` is supposed to be the *full momentum-equation residual*
(`-nu*div(grad(u)) + (u.grad)u + grad(p)`), and `grad(p)` is a term of that equation
regardless of inf-sup considerations. Cylinder's own `R1` doesn't include it either,
but legitimately not because of LBB -- Chorin's splitting genuinely excludes pressure
from the step SUPG stabilizes there (Step 1 has no pressure term in its own `F1` to
begin with). Cavity's monolithic velocity-pressure system has no such excuse: pressure
really is part of the same equation `u` is solved from.

**Fixed (2026-10-21, same day): added `grad(p)` to `R1`.** Since `R1` now depends on
both trial functions (`u` and `p`), the SUPG term no longer belongs entirely in the
velocity-velocity block -- it splits: the `u`-dependent piece
(`dot(u_n, nabla_grad(u)) - nu*div(grad(u))`) stays in `nitsche_a00` as before, and the
new `grad(p)` piece becomes its own bilinear form in `(w, p)`, added into the
velocity-pressure block alongside the existing `pressure = -p*div(w)*dx` term:
`a = [[diffusive + nitsche_a00 + convective, pressure + supg_pressure], [continuity, zero_pp]]`.

**Re-measured again (Re=400, coarse mesh) -- this fixed it**: `weak_SUPG` moved from
`-0.720` (42% low) to `-0.976` (22% low), now matching plain `weak` (`-0.969`) almost
exactly, both panels. The systematic undershoot is gone; the remaining `weak`/`weak_SUPG`
vs. `strong`/Ghia gap is just the ordinary Nitsche-penalty approximation bias (weak BCs
don't enforce the boundary condition exactly), not a SUPG artifact. Confirms the
diagnosis: a stabilization term built from a genuinely complete residual does what the
textbook guarantee promises -- `weak_SUPG` and `weak` should (and now do) agree closely,
since a correctly consistent SUPG term adds ~nothing extra at a solution that's already
smooth enough not to need stabilizing.
