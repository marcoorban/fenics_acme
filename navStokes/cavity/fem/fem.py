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
    locate_dofs_topological,
)
from dolfinx.io import VTKFile, gmsh
from dolfinx.mesh import locate_entities_boundary
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx.fem.petsc import LinearProblem
from ufl import TestFunction, TrialFunction, div, dot, dx, grad, inner, nabla_grad

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

fdim = mesh.topology.dim - 1
left = locate_entities_boundary(mesh, fdim, lambda x: np.isclose(x[0], 0))
right = locate_entities_boundary(mesh, fdim, lambda x: np.isclose(x[0], 1))
top = locate_entities_boundary(mesh, fdim, lambda x: np.isclose(x[1], 1))
bottom = locate_entities_boundary(mesh, fdim, lambda x: np.isclose(x[1], 0))

# No-slip on left/right/bottom: u = (0, 0). grad(p) = 0 there needs no
# Dirichlet BC of its own -- it's the natural (do-nothing) condition that
# already holds on any boundary with no pressure BC imposed on it.
noslip = Constant(mesh, PETSc.ScalarType((0, 0)))
lid = Constant(mesh, PETSc.ScalarType((20, 0)))
bcu = [
    dirichletbc(noslip, locate_dofs_topological(V, fdim, left), V),
    dirichletbc(noslip, locate_dofs_topological(V, fdim, right), V),
    dirichletbc(noslip, locate_dofs_topological(V, fdim, bottom), V),
    dirichletbc(lid, locate_dofs_topological(V, fdim, top), V),
]

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

nu = Constant(mesh, PETSc.ScalarType(0.01))  # Re = U*L/nu = 1*1/0.01 = 100

convective = inner(dot(u_n, nabla_grad(u_n)), w) * dx
diffusive = nu * inner(grad(u), grad(w)) * dx
pressure = -p * div(w) * dx
continuity = div(u) * q * dx

zero_pp = Constant(mesh, PETSc.ScalarType(0)) * p * q * dx  # explicit zero
# block, not None -- the pressure pin BC below needs a real (1, 1) block
# to apply to.
a = [[diffusive, pressure], [continuity, zero_pp]]
L = [-convective, Constant(mesh, PETSc.ScalarType(0)) * q * dx]

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

results_dir = Path(__file__).parent / "results" / "strong"
results_dir.mkdir(parents=True, exist_ok=True)
with VTKFile(mesh.comm, results_dir / "u.pvd", "w") as f:
    f.write_function(u_)
with VTKFile(mesh.comm, results_dir / "p.pvd", "w") as f:
    f.write_function(p_)
