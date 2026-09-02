# Navier-Stokes FEM Project — Overview

Captures the objectives, assumptions, and decisions established in the planning
conversation for this project, before any code was written. Update this file
as decisions change.

## Objective

Build a FEM solver for incompressible, viscous **3D turbulent channel flow**
(flow between two infinite parallel plates) using the Navier-Stokes
equations, implemented in FEniCS (`dolfinx`), supporting both strongly
(Dirichlet) and weakly (Nitsche) enforced boundary conditions. An OpenFOAM
simulation will serve as the baseline the FEM results are validated against.

## Governing equations

- Unsteady, incompressible Navier-Stokes — both the time-derivative and
  convective terms are retained.
- Flow regime is **3D turbulent** (DNS/LES-scale, not laminar plane
  Poiseuille), which is why those terms matter: for 2D fully-developed
  laminar channel flow the convective term is identically zero at the exact
  solution, so keeping it would add nonlinear-solver cost for no physical
  effect. In 3D turbulent channel flow both terms are physically active.
- Mean flow statistics (mean velocity profile, and any turbulence
  quantities of interest) are obtained by **time-averaging the transient
  FEM solution**, not from a RANS turbulence model.

## Boundary condition treatment

Follows the pattern already established in `2Dskew/fem.py`
(`StrongForm` / `WeakForm` / `StrongWeakForm`): every boundary-condition
entry in config can be imposed either as a `DirichletBC` or via Nitsche
terms added to the variational form, selected per-run.

## Code structure decision

- New project lives in `navier-stokes-FEM/`.
- The FEM solver itself lives in `navier-stokes-FEM/navier-stokes/`,
  seeded by copying `2Dskew/fem.py` (and, per the repo's existing
  convention of a per-directory `project_io.py`, that file too) as a
  starting scaffold — then adapted from the scalar advection-diffusion
  problem to vector incompressible Navier-Stokes.

## OpenFOAM baseline

- Does not exist yet — needs to be created from scratch.
- Uses the **openfoam.org (OpenFOAM Foundation)** line, not openfoam.com
  (OpenCFD/ESI). Currently run via a Docker container on **v13**, with a
  planned upgrade to **v14**.
- Must confirm the actual solver interface inside the container before
  writing any case dictionaries: the Foundation line modularized its
  application solvers into `foamRun` around v11, so standalone apps like
  `icoFoam` / `pimpleFoam` may not exist in v13 the way older tutorials
  assume. Check with `ls $FOAM_APPBIN` / `foamRun -listSolvers` inside the
  container first.
- Should be scoped to match the FEM side's physics (turbulent, time-
  resolved) rather than defaulting to a RANS mean-flow model, so the two
  results are comparable on the same physical basis.

## Planned build order (incremental verification ladder)

Rather than writing the full 3D turbulent Nitsche solver directly, build up
in stages, each with a known-good check before moving on:

1. Steady Stokes, strong BCs only, Taylor-Hood mixed element (2D) —
   verify against the analytic plane-Poiseuille profile
   `u(y) = (G / 2μ) · y(H − y)`.
2. Add the convective term (`NonlinearProblem` + Newton) — same analytic
   answer still applies (2D laminar), so this checks the nonlinear
   machinery in isolation.
3. Add time-stepping (still 2D laminar) — same steady answer at
   convergence, checking the transient machinery.
4. Add Nitsche (weak) boundary conditions — now with a known-good
   reference at every prior step.
5. Scale up to the real 3D turbulent case (mesh, resolution, `Re_tau`
   target, and comparison quantities still to be decided).

## Known risks / open items flagged during planning

- **Pressure datum**: if every boundary carries a velocity condition
  (strong or weak), pressure is only determined up to a constant and the
  system is singular. Needs an outlet traction / "do-nothing" condition or
  a zero-mean pressure constraint — decide when tagging boundaries, not
  when a solve fails to converge.
- **Vector Nitsche terms need the full traction, not just diffusive flux**:
  consistency and adjoint-consistency terms must use
  `σ(u, p)·n = (2ν sym(∇u) − pI)·n`, not just `-κ∇u·n` as in the scalar
  advection-diffusion case. Dropping the pressure part silently makes the
  weak BC inconsistent.
- **Inflow/outflow indicator becomes nonlinear**: in `2Dskew/fem.py`,
  `inflow()`/`outflow()` are built from a constant `beta`; here the
  advecting velocity is the unknown `u` itself, so the indicator is
  solution-dependent.
- **`params.yaml` schema change**: boundary values need to become
  vectors/expressions (e.g. a parabolic inlet profile), not the scalars
  used in the advection-diffusion case.
- **Dockerfile base image**: `2Dskew/Dockerfile` pins
  `dolfinx/dolfinx:nightly`. For a build meant to be compared against a
  fixed OpenFOAM version and reported on, a stable dolfinx tag is safer
  than a moving nightly one — not yet decided whether to change this.
- **3D turbulent mesh/compute cost is unscoped**: `Re_tau` target, mesh
  resolution (wall units), and which turbulence statistics to compare
  against OpenFOAM are not yet decided.

## Explicitly ruled out

- Treating this as 2D laminar Stokes/steady-state flow (the original
  framing implied this, but the user confirmed the goal is 3D turbulent
  channel flow with time-averaged statistics).
- OpenFOAM baseline using a RANS turbulence model while FEM stays
  laminar/time-resolved — ruled out as an apples-to-oranges comparison.
