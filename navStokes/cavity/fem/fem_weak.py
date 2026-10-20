import argparse
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
    sqrt,
    sym,
)

parser = argparse.ArgumentParser()
parser.add_argument("mesh", help="Path to the .msh file to read.")
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
U_lid = 0.005
Re = U_lid * L / nu_val  # = Pe for momentum transport (nu is the only diffusivity)

noslip = Constant(mesh, PETSc.ScalarType((0, 0)))
lid = Constant(mesh, PETSc.ScalarType((U_lid, 0)))
dirichlet_boundaries = [(1, lid), (2, noslip)]
bcu = []

# --- Steady Navier-Stokes momentum + continuity: convective (explicit,
# frozen at u_n) + diffusive + pressure-gradient + continuity. Convection
# has no dependence on the unknowns, so this is linear in (u, p) -- a
# Stokes-type system with an explicit convective load, not a Newton
# residual. Kinematic viscosity nu (no separate rho -- p is kinematic
# pressure). --- #
u = TrialFunction(V)
p = TrialFunction(Q)
w = TestFunction(V)
q = TestFunction(Q)

u_n = Function(V, name="u_n")  # previous/known iterate -- convection is
# frozen here, fully explicit, not linearized against the unknown u
u_ = Function(V, name="u")  # solution, filled in by LinearProblem.solve()
p_ = Function(Q, name="p")

nu = Constant(mesh, PETSc.ScalarType(nu_val))

convective = inner(dot(u_n, nabla_grad(u_n)), w) * dx
diffusive = nu * inner(grad(u), grad(w)) * dx
pressure = -p * div(w) * dx
continuity = div(u) * q * dx

# --- Nitsche (weak) enforcement of dirichlet_boundaries, replacing bcu,
# plus SUPG. Nitsche uses its own stress operator (2*nu*sym(grad(.))) --
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

# SUPG, stabilizing convection only: R1 = (u_n . grad)u_n. The diffusive
# (-nu*div(grad(u))) piece is deliberately left out -- feeding a P2
# field's strong-form (second-derivative) Laplacian into a bilinear form,
# scaled by an unbounded tau (no time step here to cap it, unlike the
# unsteady reference), severely ill-conditions the system and makes the
# Nitsche-enforced BC go unsatisfied (confirmed empirically: with it,
# lid mean velocity collapses to ~0.5 instead of ~20 and the solve takes
# 15x longer; without it, the solution is sane). R1 has no dependence on
# the trial function u at all now, so SUPG contributes only to L0, not
# a00 -- it's a pure stabilizing load, not a bilinear term.
a_transport = u_n  # frozen transport velocity for SUPG's tau/R1
tau = 1 / sqrt(
    (2 * sqrt(dot(a_transport, a_transport)) / h) ** 2 + (4 * nu / h**2) ** 2
)
supg_test = dot(a_transport, nabla_grad(w))
R1 = dot(u_n, nabla_grad(u_n))
nitsche_L0 += -tau * dot(supg_test, R1) * dx  # SUPG (known part)

zero_pp = Constant(mesh, PETSc.ScalarType(0)) * p * q * dx  # explicit zero
# block, not None -- the pressure pin BC below needs a real (1, 1) block
# to apply to.
a = [[diffusive + nitsche_a00, pressure], [continuity, zero_pp]]
L = [-convective + nitsche_L0, Constant(mesh, PETSc.ScalarType(0)) * q * dx]

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
# Bootstrap solve with u_n = 0 (pure Stokes, no convection), then freeze
# u_n at that result and solve once more -- this is what actually engages
# the explicit (u_n . grad) u_n term; a single step, not iterated to
# convergence, since the point here is to stay crude.
print("Starting solver...")
steps = tqdm(
    ["Stokes bootstrap (u_n = 0)", "Explicit convection (u_n = bootstrap result)"],
    desc="Solving",
)
for step in steps:
    steps.set_postfix_str(step)
    u_, p_ = problem.solve()
    u_n.x.array[:] = u_.x.array
print("Solver finished.")

results_dir = Path(__file__).parent / "results" / "weak" / f"Re_{Re:g}"
results_dir.mkdir(parents=True, exist_ok=True)
with VTKFile(mesh.comm, results_dir / "u.pvd", "w") as f:
    f.write_function(u_)
with VTKFile(mesh.comm, results_dir / "p.pvd", "w") as f:
    f.write_function(p_)
