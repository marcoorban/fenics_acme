import argparse
from pathlib import Path

import ufl
import numpy as np
from tqdm import tqdm
from dolfinx import fem
from dolfinx.io import gmsh, VTKFile
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.geometry import bb_tree, compute_collisions_points, compute_colliding_cells
from ufl import dot, grad, inner, dx, div, CellDiameter, sqrt
from mpi4py import MPI
from petsc4py import PETSc
from basix.ufl import element, mixed_element


def lower_right_corner(x):
    return np.logical_and(np.isclose(x[0], 1.0), np.isclose(x[1], 0.0))


def eval_at_points(mesh, function, points):
    """Evaluate a dolfinx Function at arbitrary (N, 3) points.

    Serial only -- assumes every point lies in a cell owned by this
    (single) process, which holds since gmshio.read_from_msh was read
    with rank=0 and no mesh partitioning is in play.
    """
    tree = bb_tree(mesh, mesh.topology.dim)
    candidates = compute_collisions_points(tree, points)
    colliding = compute_colliding_cells(mesh, candidates, points)
    cells = np.array(
        [colliding.links(i)[0] for i in range(len(points))], dtype=np.int32
    )
    return function.eval(points, cells)


parser = argparse.ArgumentParser(
    prog="Cylinder cavity solver (Newton + SUPG)",
    description="""Solves incompressible Navier-Stokes equation in cavity flow
            using FEM method, via Newton's method (SNES) on the full
            nonlinear residual, with SUPG convection stabilization added
            (no PSPG -- Taylor-Hood P2/P1 already satisfies LBB).""",
    usage="fem_SUPG.py mesh Re results_dir",
)
parser.add_argument("mesh", help="Path to the .msh file to read.")
parser.add_argument(
    "reynolds", type=float, help="Reynolds number to use in simulation."
)
parser.add_argument("results_dir", help="Directory to write results into.")
args = parser.parse_args()

assert args.reynolds >= 1.0

TARGET_RE = args.reynolds
# Matches jupyter/Steady_Lid_Cavity_SUPG_PSPG_FEniCSx_20261022.ipynb's own
# continuation schedule -- RE_STEPS = [1, 10, 50, 100, 200, TARGET_RE] for
# its TARGET_RE=400 case -- generalized to any target: those fixed low-Re
# anchors, then doubling (not fem_newton.py's tripling) each step after,
# so every SUPG nonlinear solve starts closer to the previous stage.
RE_STEPS = [re for re in (1.0, 10.0, 50.0, 100.0, 200.0) if re < TARGET_RE] or [
    TARGET_RE
]
while RE_STEPS[-1] < TARGET_RE:
    RE_STEPS.append(min(RE_STEPS[-1] * 2.0, TARGET_RE))

# Ensure that the last step is the one we want and that
# they all go in ascending order
assert RE_STEPS[-1] == TARGET_RE
assert all(b > a for a, b in zip(RE_STEPS, RE_STEPS[1:]))

# Read the mesh
meshData = gmsh.read_from_msh(args.mesh, MPI.COMM_WORLD, rank=0, gdim=2)
mesh = meshData.mesh
ft = meshData.facet_tags

# velocity element space is one order higher than P to preserve
# inf-sup (LBB) condition.
v_cg = element("Lagrange", mesh.basix_cell(), 2, shape=(mesh.geometry.dim,))
p_cg = element("Lagrange", mesh.basix_cell(), 1)
# V = fem.functionspace(mesh, v_cg)
# Q = fem.functionspace(mesh, p_cg)

W = fem.functionspace(mesh, mixed_element([v_cg, p_cg]))
V, _ = W.sub(0).collapse()
Q, _ = W.sub(1).collapse()
W0 = W.sub(0)
W1 = W.sub(1)

# Single mixed unknown Function -- Newton solves F(solution) = 0
# directly against the true nonlinear residual, so (unlike the Picard
# solver in fem_strong_relaxed.py) there's no separate trial function
# or lagged u_n: solution IS the current Newton iterate throughout, and
# SNES's own loop (not a hand-rolled Python one) drives it to the root.
solution = fem.Function(W)
solution.name = "velocity_pressure"
u, p = ufl.split(solution)
v, q = ufl.TestFunctions(W)

# Physical constants
length = 1.0
nuVal = 1.0
nu = fem.Constant(mesh, PETSc.ScalarType(1.0))

# Boundary conditions
fdim = mesh.topology.dim - 1

# Lid velocity
# A vector-valued Dirichlet BC on a sub-space of a genuine mixed_element
# space can't use a plain Constant here: W0's reported block size is 1
# (dolfinx 0.11), which doesn't match a size-2 Constant. dolfinx's own
# fix for this is to use a Function living in the collapsed space V
# instead, paired with the (W0, V) dof mapping.
lid_fn = fem.Function(V)
lid_facets = ft.find(
    1
)  # This 1 is a magic number! It is the first tag of the gmsh mesh, which has been assigned
# to the lid.
dofs_lid = fem.locate_dofs_topological((W0, V), fdim, lid_facets)
bc_lid = fem.dirichletbc(lid_fn, dofs_lid, W0)

# NoSlip boundary condition -- stays zero for the whole run, so the
# default zero-initialized Function needs no further updates.
noSlip_fn = fem.Function(V)
wall_facets = ft.find(2)
dofs_wall = fem.locate_dofs_topological((W0, V), fdim, wall_facets)
bc_wall = fem.dirichletbc(noSlip_fn, dofs_wall, W0)

# Pressure boundary condition
# A subspace of a mixed_element space can't tabulate its own
# coordinates, so locate_dofs_geometrical needs the (W1, Q) paired form
# -- but a Constant-valued BC needs a single dof array (see the lid/wall
# BCs above), so take just the W1-numbered column back out of the pair.
zeroPressure = fem.Constant(mesh, PETSc.ScalarType(0))
dofs_corner_pair = fem.locate_dofs_geometrical((W1, Q), lower_right_corner)
dofs_corner = dofs_corner_pair[0]
bc_pressure = fem.dirichletbc(zeroPressure, dofs_corner, W1)

bcs = [bc_lid, bc_wall, bc_pressure]

# Forcing term -- it's just zero but adding here for completeness
f = fem.Constant(mesh, PETSc.ScalarType((0, 0)))

# Navier-Stokes terms, strongly enforced boundary conditions. Same weak
# form as fem_strong_relaxed.py, except convection now uses the real
# unknown u on both sides (genuinely nonlinear/quadratic in u) instead
# of a lagged u_n -- Newton solves this directly rather than
# linearizing it by hand.
#
# dot(grad(u), u), NOT dot(u, grad(u)): UFL's grad(u)[i,j] = d(u_i)/
# d(x_j), so (u.grad)u (component j = sum_i u_i * d(u_j)/d(x_i)) is the
# matrix-vector contraction dot(grad(u), u), contracting grad(u)'s
# SECOND index against u. dot(u, grad(u)) contracts the FIRST index
# instead, giving sum_i u_i * d(u_i)/d(x_j) = grad(0.5*dot(u,u)) -- a
# pure gradient field, not a transport term. Verified numerically
# against a hand-computed (u.grad)u on an exactly-representable linear
# field. Since a pure-gradient "convection" term is fully absorbable
# into the pressure unknown (only grad(p) appears, p is otherwise
# unconstrained), that bug made every Re collapse onto the Re=1 Stokes
# velocity field exactly rescaled -- confirmed empirically: the
# Re=1 and Re=400 u_vertical.dat profiles were identical once divided
# by their respective lid speeds, to 5 decimal places.
convection = inner(dot(grad(u), u), v) * dx
pressure = inner(-grad(p), v) * dx
diffusion = nu * inner(grad(u), grad(v)) * dx
continuity = -inner(q, div(u)) * dx
force = inner(f, v) * dx
F = convection + diffusion - pressure + continuity - force

# SUPG convection stabilization. Unlike fem_newton.py, this is meant for
# the regime where the mesh can't resolve convection relative to
# diffusion (high cell Peclet number, see plot_peclet.py) -- confirmed
# to actually matter on a deliberately under-resolved 10x10 structured
# mesh at Re=1000, which still converged in Newton's sense but produced
# a spurious overshoot near the lid and a wildly unphysical pressure
# field. No PSPG term: that one exists to patch equal-order (e.g. P1/P1)
# elements violating LBB, and we're already on Taylor-Hood (P2/P1),
# which doesn't need it.
#
# Strong-form momentum residual -- zero for the exact solution, nonzero
# for any discrete approximation. This is what SUPG tests against, which
# is what makes it a *consistent* stabilization (it vanishes exactly at
# the true solution, so it doesn't change what equation is being
# solved, only how well-behaved the discrete approximation is).
velocity_laplacian = div(grad(u))
R = dot(grad(u), u) - nu * velocity_laplacian + grad(p) - f

# CellDiameter is the max distance between any two vertices in the cell
# -- unlike a hand-derived triangle-specific formula (e.g. from
# jupyter/Steady_Lid_Cavity_SUPG_PSPG_FEniCSx_20261022.ipynb, which
# assumes triangles), this is exact for any convex cell shape, including
# the quadrilaterals this mesh actually uses.
h = CellDiameter(mesh)
# sqrt(dot(u,u))'s derivative is u/|u|, singular at u=0 -- and the
# Newton initial guess starts at solution=0 everywhere, so the very
# first Jacobian evaluation would divide by zero without this floor.
# Negligible next to any genuine velocity scale here (O(1)-O(1e4)).
speed = sqrt(dot(u, u) + 1.0e-12)
# Tezduyar/Franca-style blend: tau ~ h/(2|u|) when convection dominates
# (high cell Peclet number), tau ~ h^2/(4*nu) when diffusion dominates --
# robust across both regimes without needing a coth()-based 1D-optimal
# formula or any triangle-specific geometry.
tau_supg = 1.0 / sqrt((2.0 * speed / h) ** 2 + (4.0 * nu / h**2) ** 2)

# Test function perturbed by tau*(u.grad)v, same advective-derivative
# convention as the convection term above (grad(v) first, contracting
# its second index against u).
supg = tau_supg * inner(dot(grad(v), u), R) * dx
F = F + supg

# Newton solver settings, matching
# jupyter/Steady_Lid_Cavity_SUPG_PSPG_FEniCSx_20261022.ipynb's
# proven-fast setup: backtracking line search for globalization, and a
# direct LU factorization for each Newton step's linear solve (no
# iterative-KSP residual noise contaminating convergence, unlike
# LinearProblem's PETSc defaults in the Picard version).
petsc_options = {
    "snes_type": "newtonls",
    "snes_linesearch_type": "bt",
    "snes_atol": 1.0e-10,
    "snes_rtol": 1.0e-8,
    "snes_max_it": 40,
    "ksp_type": "preonly",
    "pc_type": "lu",
}

problem = NonlinearProblem(
    F,
    solution,
    bcs=bcs,
    petsc_options_prefix="cavity_newton_",
    petsc_options=petsc_options,
)


def newton_monitor(snes, iteration, residual_norm):
    tqdm.write(f"  Newton {iteration:2d}: |F| = {residual_norm:.6e}")


problem.solver.setMonitor(newton_monitor)

re_progress = tqdm(RE_STEPS, desc="Re continuation", unit="step")
for re in re_progress:
    re_progress.set_postfix(Re=re)
    tqdm.write(f"Solving at Re = {re:g}")
    U_lid = re * nuVal / length
    lid_fn.interpolate(
        lambda x: np.vstack([np.full(x.shape[1], U_lid), np.zeros(x.shape[1])])
    )
    problem.solve()  # reuses the previous Re stage's solution as the initial guess
    reason = problem.solver.getConvergedReason()
    iterations = problem.solver.getIterationNumber()
    tqdm.write(f"  Newton finished: {iterations} iteration(s), reason={reason}")
    if reason <= 0:
        raise RuntimeError(
            f"Newton solve did not converge at Re={re:g} (SNES reason {reason})."
        )

solution.x.scatter_forward()
uh = solution.sub(0).collapse()

# Write results
ph = solution.sub(1).collapse()
uh.name = "u"
ph.name = "p"

results_dir = Path(args.results_dir)
dat_dir = results_dir / "dat"
dat_dir.mkdir(parents=True, exist_ok=True)

# ParaView-readable output.
with VTKFile(mesh.comm, results_dir / "u.pvd", "w") as vtk_file:
    vtk_file.write_function(uh)
with VTKFile(mesh.comm, results_dir / "p.pvd", "w") as vtk_file:
    vtk_file.write_function(ph)

# Midline profiles, matching Ghia, Ghia & Shin (1982)'s own convention:
# u along the VERTICAL centerline (x=0.5), v along the HORIZONTAL
# centerline (y=0.5) -- see navStokes/cavity/benchmarks/.
n_sample = 129
y_line = np.linspace(0.0, 1.0, n_sample)
x_line = np.linspace(0.0, 1.0, n_sample)

vertical_points = np.column_stack([0.5 * np.ones(n_sample), y_line, np.zeros(n_sample)])
horizontal_points = np.column_stack(
    [x_line, 0.5 * np.ones(n_sample), np.zeros(n_sample)]
)

u_vertical = eval_at_points(mesh, uh, vertical_points)[:, 0]
v_horizontal = eval_at_points(mesh, uh, horizontal_points)[:, 1]

np.savetxt(
    dat_dir / "u_vertical.dat",
    np.column_stack([y_line, u_vertical]),
    header="y u",
    comments="# ",
)
np.savetxt(
    dat_dir / "v_horizontal.dat",
    np.column_stack([x_line, v_horizontal]),
    header="x v",
    comments="# ",
)

# Full velocity + pressure field on a regular grid, for streamline/
# contour plotting with plain matplotlib (no pyvista/vtk/scipy needed --
# none of those besides vtk itself are available in this environment,
# and vtk's API is awkward for this compared to just sampling the
# Functions we already have in memory).
n_grid = 101
xs = np.linspace(0.0, 1.0, n_grid)
ys = np.linspace(0.0, 1.0, n_grid)
Xg, Yg = np.meshgrid(xs, ys)
grid_points = np.column_stack([Xg.ravel(), Yg.ravel(), np.zeros(Xg.size)])

uv_grid = eval_at_points(mesh, uh, grid_points)
p_grid = eval_at_points(mesh, ph, grid_points)

np.savez(
    results_dir / "field_grid.npz",
    x=xs,
    y=ys,
    u=uv_grid[:, 0].reshape(Xg.shape),
    v=uv_grid[:, 1].reshape(Xg.shape),
    p=p_grid[:, 0].reshape(Xg.shape),
    Re=TARGET_RE,
)

print(f"Results written to {results_dir}")
