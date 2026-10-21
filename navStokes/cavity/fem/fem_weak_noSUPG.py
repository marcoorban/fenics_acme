import argparse
import time
from pathlib import Path

import numpy as np
from tqdm import tqdm
from basix.ufl import element
from dolfinx.fem import (
    Constant,
    Function,
    functionspace,
    dirichletbc,
    locate_dofs_geometrical,
)
from dolfinx.io import VTKFile, gmsh
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx.fem.petsc import LinearProblem
from ufl import (
    CellDiameter,
    FacetNormal,
    Measure,
    TestFunction,
    TrialFunction,
    div,
    dot,
    dx,
    grad,
    inner,
    nabla_grad,
    sym,
)

parser = argparse.ArgumentParser()
parser.add_argument("mesh", help="Path to the .msh file to read.")
parser.add_argument(
    "reynolds", type=float, help="Reynolds number Re = U_lid*L/nu (L=1, nu=0.01)."
)
args = parser.parse_args()

meshData = gmsh.read_from_msh(args.mesh, MPI.COMM_WORLD, rank=0, gdim=2)
mesh = meshData.mesh
ft = meshData.facet_tags

v_cg = element("Lagrange", mesh.basix_cell(), 2, shape=(mesh.geometry.dim,))
s_cg = element("Lagrange", mesh.basix_cell(), 1)
V = functionspace(mesh, v_cg)
Q = functionspace(mesh, s_cg)

# Velocity BCs are enforced weakly (Nitsche) below instead of strongly, so
# bcu is empty -- these (tag, target) pairs are consumed by the Nitsche
# terms via ds measures on ft, not by dirichletbc()/locate_dofs_topological.
# Tags match mesh.py's physical groups: lid=1 (moving lid), wall=2 (the
# other three, fixed, no-slip edges).
L = 1.0  # unit square cavity
nu_val = 0.01
Re = args.reynolds  # = Pe for momentum transport (nu is the only diffusivity)
U_lid = Re * nu_val / L
lid_per_re = nu_val / L  # captured now -- L gets reassigned to the RHS
# block list further down, so this ratio can't be recomputed from L later

noslip = Constant(mesh, PETSc.ScalarType((0, 0)))
lid = Constant(mesh, PETSc.ScalarType((U_lid, 0)))
dirichlet_boundaries = [(1, lid), (2, noslip)]
bcu = []

# --- Steady Navier-Stokes momentum + continuity: convective + diffusive
# + pressure-gradient + continuity. The convective term (u.grad)u is
# nonlinear (it is quadratic in the unknown u), so it is Picard-
# linearized: one copy of u is replaced by the known/lagged u_n, exactly
# like a passive advection-diffusion problem where u_n plays the role of
# the prescribed transport velocity "a" and the remaining trial function
# u plays the role of the transported quantity "f". This makes each
# individual solve genuinely linear in (u, p); repeatedly solving and
# setting u_n <- u (see the Picard loop below) drives u_n -> u, at which
# point the Picard-linearized equation coincides with the true nonlinear
# steady Navier-Stokes equation. Kinematic viscosity nu (no separate rho
# -- p is kinematic pressure). --- #
u = TrialFunction(V)
p = TrialFunction(Q)
w = TestFunction(V)
q = TestFunction(Q)

u_n = Function(V, name="u_n")  # known/lagged transport velocity -- updated
# after every Picard iteration, not just once
u_ = Function(V, name="u")  # solution, filled in by LinearProblem.solve()
p_ = Function(Q, name="p")

nu = Constant(mesh, PETSc.ScalarType(nu_val))

# Picard linearization: u_n (known) transports the real trial function u
# (unknown), so this term is bilinear in (u, w) for a given u_n and
# belongs in the system matrix, not the right-hand side.
convective = inner(dot(u_n, nabla_grad(u)), w) * dx
diffusive = nu * inner(grad(u), grad(w)) * dx
pressure = -p * div(w) * dx
continuity = div(u) * q * dx

# --- Nitsche (weak) enforcement of dirichlet_boundaries, replacing bcu.
# No SUPG here (see fem_weak.py for that variant) -- this is otherwise
# identical. Nitsche uses its own stress operator (2*nu*sym(grad(.))) --
# a self-contained, standard Nitsche formulation, independent of the bulk
# diffusive term above (which stays a plain Laplacian, as in fem.py).
# Derived from the residual form ibp(u,w) + consistency(u,w,g) +
# penalty(u,w,g), with consistency/penalty expanded in (u - g) and split
# into their bilinear-in-u piece (-> a00) and known-g piece (-> L0, with
# the sign flip that moving a residual term to the RHS of a*x=L implies).
# F1/lhs()/rhs() don't work here -- V and Q are separate function spaces,
# not one MixedFunctionSpace, so a combined residual has two different
# Arguments at the same number/part and UFL can't split it; the a00/L0
# block below is built by hand instead.
# --- #
vDegree = 2  # V is P2


def sigma_u(uu):
    return 2 * nu * sym(grad(uu))


n = FacetNormal(mesh)
h = CellDiameter(mesh)
# Dimensionless Nitsche penalty constant; the actual penalty coefficient
# is beta*nu/h (not beta/h) -- it must scale with nu to match the ibp/
# consistency terms, which carry nu via sigma_u. Omitting nu here made
# the penalty ~100x too stiff relative to the viscous terms (nu=0.01),
# swamping the system and degrading the solution.
beta = Constant(mesh, PETSc.ScalarType(10 * vDegree**2))  # tune if unstable

nitsche_a00 = 0
nitsche_L0 = 0
for tag, g in dirichlet_boundaries:
    ds_tag = Measure("ds", domain=mesh, subdomain_data=ft, subdomain_id=tag)
    nitsche_a00 += -dot(w, dot(sigma_u(u), n)) * ds_tag  # ibp
    nitsche_a00 += -dot(dot(sigma_u(w), n), u) * ds_tag  # consistency (bilinear part)
    nitsche_L0 += -dot(dot(sigma_u(w), n), g) * ds_tag  # consistency (known part)
    nitsche_a00 += (beta * nu / h) * dot(u, w) * ds_tag  # penalty (bilinear part)
    nitsche_L0 += (beta * nu / h) * dot(g, w) * ds_tag  # penalty (known part)

zero_pp = Constant(mesh, PETSc.ScalarType(0)) * p * q * dx  # explicit zero
# block, not None -- the pressure pin BC below needs a real (1, 1) block
# to apply to.
a = [[diffusive + nitsche_a00 + convective, pressure], [continuity, zero_pp]]
L = [nitsche_L0, Constant(mesh, PETSc.ScalarType(0)) * q * dx]

# Pressure has no Dirichlet BC anywhere else (grad(p)=0 is natural on
# every wall here), so it's only defined up to an additive constant --
# the saddle-point matrix is singular by one rank. Pin it at a single
# point (the (0,0) corner) to remove that nullspace.
bcp = dirichletbc(
    PETSc.ScalarType(0),
    locate_dofs_geometrical(Q, lambda x: np.isclose(x[0], 0) & np.isclose(x[1], 0)),
    Q,
)

problem = LinearProblem(
    a, L, bcs=bcu + [bcp], u=[u_, p_], petsc_options_prefix="cavity_"
)

picard_tol = 5e-6
picard_max_iter = 200
picard_omega = 0.5

# Reynolds-number continuation: ramp from a moderate Re up to the target
# Re through a geometric sequence (doubling each stage), reusing the
# previous stage's converged u_n as the next stage's starting point
# instead of jumping straight from rest to the target Re. Plain Picard
# starting from rest stalls in a low-amplitude oscillation at Re=400
# instead of converging (confirmed empirically: 200 iterations,
# rel_change plateaued around 2e-5, never reached picard_tol); starting
# each stage much closer to its own fixed point avoids that. Easy cases
# (Re<=50) reduce to a single direct stage, so this adds no overhead
# there. See README.md for the full scheme/parameter writeup.
continuation_steps = []
re_stage = min(50.0, Re)
while re_stage < Re:
    continuation_steps.append(re_stage)
    re_stage = min(re_stage * 2, Re)
continuation_steps.append(Re)

print("Starting solver...")
t_start = time.time()
total_iterations = 0
for stage_re in continuation_steps:
    lid.value = np.array([stage_re * lid_per_re, 0.0], dtype=PETSc.ScalarType)
    progress = tqdm(range(1, picard_max_iter + 1), desc=f"Picard (Re={stage_re:g})")
    for it in progress:
        u_, p_ = problem.solve()
        rel_change = float(
            np.linalg.norm(u_.x.array - u_n.x.array) / np.linalg.norm(u_.x.array)
        )
        u_n.x.array[:] = u_n.x.array + picard_omega * (u_.x.array - u_n.x.array)
        progress.set_postfix_str(f"rel_change={rel_change:.2e}")
        if rel_change < picard_tol:
            break
    progress.close()
    total_iterations += it
    if rel_change >= picard_tol:
        print(f"Picard did NOT converge at Re={stage_re:g} in {picard_max_iter} "
              f"iterations (rel_change={rel_change:.2e})")
t_elapsed = time.time() - t_start
print(f"Solver finished after {total_iterations} total Picard iteration(s) "
      f"across {len(continuation_steps)} continuation stage(s), "
      f"rel_change={rel_change:.2e}, time={t_elapsed:.1f}s")

mesh_stem = Path(args.mesh).stem
results_dir = Path(__file__).parent / "results" / "weak_noSUPG" / mesh_stem / f"Re_{Re:g}"
results_dir.mkdir(parents=True, exist_ok=True)
with VTKFile(mesh.comm, results_dir / "u.pvd", "w") as f:
    f.write_function(u_)
with VTKFile(mesh.comm, results_dir / "p.pvd", "w") as f:
    f.write_function(p_)

# Append this run's performance to the shared cross-case log (one line
# per solver invocation, across every method/mesh/Re combination).
timing_path = Path(__file__).parent / "results" / "dat" / "solver_timing.dat"
timing_path.parent.mkdir(parents=True, exist_ok=True)
write_header = not timing_path.exists()
with open(timing_path, "a") as f:
    if write_header:
        f.write("# mesh method Re iterations time_s omega converged\n")
    f.write(
        f"{mesh_stem} weak_noSUPG {Re:g} {total_iterations} {t_elapsed:.3f} "
        f"{picard_omega:g} {int(rel_change < picard_tol)}\n"
    )
