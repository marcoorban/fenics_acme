# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A FEniCSx (`dolfinx`) finite-element solver for the 2D skew advection-diffusion
equation, stabilized with SUPG, on a rectangular domain whose left edge is
split into two boundary tags. It is the modular refactor of the standalone
prototype scripts in the sibling directories `../allStrong`, `../allWeak`,
and `../strongWeak` (each of those is a flat, single-file script exploring
one boundary-condition treatment; this directory consolidates all three into
one class hierarchy driven by config).

`nitsche_strong.py` in this directory is the original flat prototype script
kept for reference — it is not imported by `fem.py` and is not part of the
modular solver.

## Environment

The project runs against a conda environment named `fenicsx-env`
(`/opt/miniconda3/envs/fenicsx-env`, see `../../pyrightconfig.json`), which
provides `dolfinx`, `ufl`, `mpi4py`, `petsc4py`, `gmsh`, and `yaml`. Activate
it before running any script here.

## Running

```bash
python fem.py
```

This reads `params.yaml` and `rectangle.msh` (both hardcoded in
`fem.py`'s `__main__`), solves the problem with all three boundary
treatments (`StrongForm`, `WeakForm`, `StrongWeakForm`), and writes
`solutionStrong.xdmf`, `solutionWeak.xdmf`, `solutionStrongWeak.xdmf`, each
printing the min/max of the solution field.

To regenerate the mesh (`rectangle.msh`), a structured quad grid of the unit
rectangle with the left edge split at `y_split`:

```bash
python structured_mesh.py [-o out.msh] [--nx NX] [--ny NY] [--lx LX] [--ly LY] [--y-split YS]
```

Defaults (`nx=ny=65`, `y_split=0.4`) differ from `mesh.yaml` (`nx=ny=20`),
which is not currently read by any script — `structured_mesh.py` is invoked
via its own CLI flags, not `mesh.yaml`.

There is no test suite, linter, or build step in this repository.

## Architecture

`fem.py` centers on `FEM_Solver`, a base class whose `__init__` runs a fixed
pipeline: read params → read mesh → build geometry (`dx`/`ds` measures,
facet normal, cell diameter) → build the function space → build `Physics` →
build the SUPG stabilization (`tau`, `supgL_w`). Three subclasses each
override `applyStrongBCs()` / `applyWeakBCs()` to select how boundary
conditions from `params.yaml` are imposed:

- **`StrongForm`** — every entry in `boundary_conditions.strong` +
  `boundary_conditions.weak` becomes a `DirichletBC` (the split is ignored;
  see `allBCs()`).
- **`WeakForm`** — every entry from both lists is imposed weakly via
  Nitsche terms added to the bilinear/linear forms.
- **`StrongWeakForm`** — respects the YAML split: `boundary_conditions.strong`
  entries become `DirichletBC`s, `boundary_conditions.weak` entries get
  Nitsche terms.

Boundary conditions are not tagged inlet/outlet in `params.yaml`; instead,
`inflow()`/`outflow()` classify each boundary point per-quadrature-point from
the sign of `dot(beta, n)`, so a single tagged edge can be partly inlet and
partly outlet (see the docstring on `nitscheTerms()` for the upwinding logic
this drives).

`Physics` (also in `fem.py`) holds the raw YAML values for `kappa`, `alpha`,
`beta_mag`, `f` alongside the `dolfinx.fem.Constant`s built from them, so a
solver only needs one attribute (`self.physics`) instead of one per
coefficient.

Config is split across two YAML files: `params.yaml` (physics, boundary
conditions, Nitsche parameters `gamma`/`C`, polynomial family/order) loaded
via `project_io.readParams`, and `mesh.yaml` (grid resolution for
`structured_mesh.py`, currently unused by `fem.py`).

To add a new boundary-condition strategy, subclass `FEM_Solver` and override
`applyStrongBCs()`/`applyWeakBCs()` — the bilinear/linear form assembly,
SUPG term, and Nitsche term construction (`nitscheTerms()`) are shared and
should not need duplicating.
