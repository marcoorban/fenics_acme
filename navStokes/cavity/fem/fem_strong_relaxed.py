import argparse
from pathlib import Path

import ufl
import numpy as np
from tqdm import tqdm
from dolfinx import fem
from dolfinx.io import gmsh, VTKFile
from dolfinx.fem.petsc import LinearProblem
from dolfinx.geometry import bb_tree, compute_collisions_points, compute_colliding_cells
from ufl import dot, grad, inner, dx, div
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
    prog="Cylinder cavity solver",
    description="""Solves incompressible Navier-Stokes equation in cavity flow
            using FEM method""",
    usage="fem.py mesh Re",
)
parser.add_argument("mesh", help="Path to the .msh file to read.")
parser.add_argument(
    "reynolds", type=float, help="Reynolds number to use in simulation."
)
args = parser.parse_args()

assert args.reynolds >= 1.0

TARGET_RE = args.reynolds
RE_STEPS = [1.0]
while RE_STEPS[-1] < TARGET_RE:
    RE_STEPS.append(min(RE_STEPS[-1] * 3, TARGET_RE))

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

# Unknown functions that PETSc solves for.
u, p = ufl.TrialFunctions(W)
v, q = ufl.TestFunctions(W)

# Placeholder for the solution of u of the previous iteration.
u_n = fem.Function(V)
# Function that receives each solve's results
wh = fem.Function(W)

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

# Navier-Stokes terms, strongly enforced boundary conditions. Convection
# is Picard-linearized (u_n lagged, see the under-relaxed update below).
# LHS Navier-Stokes
convection = inner(dot(u_n, grad(u)), v) * dx
# RHS Navier-Stokes
pressure = inner(-grad(p), v) * dx
diffusion = nu * inner(grad(u), grad(v)) * dx
# Continuity equation
continuity = -inner(q, div(u)) * dx
# Forcing term
force = inner(f, v) * dx
a = convection + diffusion - pressure + continuity
L = force

# Solver settings
problem = LinearProblem(a, L, bcs=bcs, u=wh, petsc_options_prefix="cavity_")

max_picard_iters = 50
picard_tol = 0.1 / 100
# Under-relaxation factor on the u_n update. Plain (omega=1.0) Picard
# with lagged convection oscillates in a low-amplitude limit cycle
# instead of converging -- confirmed empirically here (Re=1 plateaued
# around 1e-3, Re=100 bounced as high as 0.3, neither reaching
# picard_tol in 50 iterations). Blending the update damps that
# oscillation into a genuine contraction; lower this (e.g. 0.5, then
# 0.3) if a given Re doesn't converge within max_picard_iters.
picard_omega = 0.5

re_progress = tqdm(RE_STEPS, desc="Re continuation", unit="step")
for re in re_progress:
    re_progress.set_postfix(Re=re)
    tqdm.write(f"Solving at Re = {re:g}")
    U_lid = re * nuVal / length
    lid_fn.interpolate(
        lambda x: np.vstack([np.full(x.shape[1], U_lid), np.zeros(x.shape[1])])
    )
    for picard_iter in range(max_picard_iters):
        problem.solve()
        uh = wh.sub(0).collapse()  # Velocity component of the new iterate
        # convergence check -- measures the raw (unrelaxed) Picard step
        # size, i.e. how far the fixed-point iteration still has to go,
        # independent of how much of that step we actually take below.
        diff_form = fem.form(inner(uh - u_n, uh - u_n) * dx)
        norm_form = fem.form(inner(uh, uh) * dx)
        diff_norm = np.sqrt(fem.assemble_scalar(diff_form))
        sol_norm = np.sqrt(fem.assemble_scalar(norm_form))
        rel_diff = diff_norm / sol_norm if sol_norm > 0 else diff_norm
        converged = rel_diff < picard_tol
        # Relaxed update: blend toward uh instead of jumping straight to
        # it, instead of u_n.x.array[:] = uh.x.array (omega=1.0).
        u_n.x.array[:] = u_n.x.array + picard_omega * (uh.x.array - u_n.x.array)
        # Printed unconditionally (not just on failure) so the last
        # value -- whichever iteration breaks the loop below, converged
        # or not -- is always visible before moving on to the next Re.
        status = "converged" if converged else ""
        tqdm.write(
            f"  iter {picard_iter + 1}/{max_picard_iters}: "
            f"relative difference = {rel_diff:.3%} ({status})"
        )
        if converged:
            break

# Write results
ph = wh.sub(1).collapse()
uh.name = "u"
ph.name = "p"

mesh_label = Path(args.mesh).stem
results_dir = Path(__file__).parent / "results" / f"Re_{TARGET_RE:g}" / mesh_label
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
