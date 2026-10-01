import os
import numpy as np
import matplotlib.pyplot as plt
import tqdm.autonotebook
import argparse
from pathlib import Path
from mpi4py import MPI
from petsc4py import PETSc

from basix.ufl import element
from dolfinx.io import XDMFFile, gmsh, VTXWriter

from dolfinx.fem import (
    Constant,
    Function,
    assemble,
    functionspace,
    assemble_scalar,
    dirichletbc,
    extract_function_spaces,
    form,
    locate_dofs_topological,
)

from dolfinx.fem.petsc import (
    apply_lifting,
    assemble_matrix,
    assemble_vector,
    create_vector,
    create_matrix,
    set_bc,
)

from ufl import (
    FacetNormal,
    Measure,
    TestFunction,
    TrialFunction,
    as_vector,
    div,
    dot,
    dx,
    inner,
    lhs,
    grad,
    nabla_grad,
    rhs,
)

from dolfinx.geometry import bb_tree, compute_collisions_points, compute_colliding_cells

# DFG benchmark case definitions -- single source of truth for each case's
# inlet peak velocity, whether it's time-ramped, and how long to run. Ubar
# (used below in the drag/lift normalization constant) is derived from Um
# here rather than hardcoded separately, so it can't drift out of sync with
# whichever case is actually selected -- that desync (Ubar implicitly fixed
# at 1.0 regardless of case) was the bug this replaces: case 1's Cd/Cl came
# out ~25x too small silently, no crash. T is case-dependent too: case 1
# only needs to run long enough to reach steady state (OpenFOAM's
# residualControl converged it at t~2.36s; 3.0s gives headroom without
# marching needlessly far past convergence -- was 8.0 unconditionally, i.e.
# >4x longer than necessary), case 2 matches the OpenFOAM 2D-2 endTime=5s
# convention, case 3 needs the full 8s to capture the whole ramp-up/down.
CASES = {
    1: {"name": "2D-1_steady", "Um": 0.3, "ramped": False, "T": 3.0},
    2: {"name": "2D-2_unsteady", "Um": 1.5, "ramped": False, "T": 5.0},
    3: {"name": "2D-3_unsteady_ramp", "Um": 1.5, "ramped": True, "T": 8.0},
}

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--inletBC",
    type=int,
    default=3,
    choices=sorted(CASES),
    help="; ".join(f"{k}. {v['name']}" for k, v in CASES.items()),
)
parser.add_argument(
    "-m",
    "--mesh",
    type=str,
    default="quad2D_2504cells.msh",
    help="Mesh file to solve on -- see cylinder/mesh/meshes_summary.csv for "
    "the catalog of available meshes and their properties. A bare filename "
    "(e.g. 'prism3D_62373cells.msh') resolves against cylinder/mesh/; an "
    "absolute or ./-relative path is used as-is. Results are kept in their "
    "own per-mesh subfolder (named after the mesh's filename stem) so "
    "different meshes don't overwrite each other's results.",
)
args = parser.parse_args()
case = CASES[args.inletBC]
Um = case["Um"]
Ubar = (2 / 3) * Um  # mean velocity of a parabolic (Poiseuille) inlet profile

# cylinder/ was reorganized 2026-10-01: meshes now live in a shared
# cylinder/mesh/ (per-method mesh/results split), not alongside each
# method's solver script.
cylinderDir = Path(__file__).parent.parent.parent

# A bare filename resolves against cylinder/mesh/ (the common case: picking
# between meshes that already live there, e.g. -m finemsh.msh); a path
# containing a directory component (absolute, or explicitly "./foo.msh") is
# used as-is instead, rather than silently resolving against the *launch*
# cwd -- that exact relative-path ambiguity bit this script once already
# (a stale duplicate mesh one level up got picked up silently).
meshArg = Path(args.mesh)
meshFile = meshArg if meshArg.parent != Path(".") else cylinderDir / "mesh" / meshArg
meshData = gmsh.read_from_msh(meshFile, MPI.COMM_WORLD, rank=0, gdim=2)
mesh = meshData.mesh
assert meshData.facet_tags is not None
ft = meshData.facet_tags
ft.name = "Facet markers"


# --- Physical and discretization parameters ---#
t = 0.0
T = case["T"]  # Final time -- case-dependent, see CASES above
dt = 1 / 1600
num_steps = int(T / dt)
k = Constant(mesh, PETSc.ScalarType(dt))
mu = Constant(mesh, PETSc.ScalarType(0.001))  # Dynamics viscosity
rho = Constant(mesh, PETSc.ScalarType(1))
D = 0.1  # cylinder diameter, matches cylinder/mesh/quad2D.py's D

# --- Boundary conditions --- #
v_cg2 = element("Lagrange", mesh.basix_cell(), 2, shape=(mesh.geometry.dim,))
s_cg1 = element("Lagrange", mesh.basix_cell(), 1)
V = functionspace(mesh, v_cg2)
Q = functionspace(mesh, s_cg1)

fdim = mesh.topology.dim - 1
gdim = 2

boundary_mapping = {"inlet": 1, "outlet": 2, "walls": 3, "cylinder": 4}


class InletVelocity:
    def __init__(self, t, Um, ramped):
        self.t = t
        self.Um = Um
        self.ramped = ramped

    def __call__(self, x):
        values = np.zeros((gdim, x.shape[1]), dtype=PETSc.ScalarType)
        peak = self.Um * np.sin(self.t * np.pi / 8) if self.ramped else self.Um
        values[0] = 4 * peak * x[1] * (0.41 - x[1]) / (0.41**2)
        return values


# Inlet
u_inlet = Function(V)
inlet_velocity = InletVelocity(t, Um, case["ramped"])
u_inlet.interpolate(inlet_velocity)
bcu_inflow = dirichletbc(
    u_inlet, locate_dofs_topological(V, fdim, ft.find(boundary_mapping["inlet"]))
)

# Walls
u_nonslip = np.array((0,) * mesh.geometry.dim, dtype=PETSc.ScalarType)
bcu_walls = dirichletbc(
    u_nonslip, locate_dofs_topological(V, fdim, ft.find(boundary_mapping["walls"])), V
)

# Cylinder
bcu_cylinder = dirichletbc(
    u_nonslip,
    locate_dofs_topological(V, fdim, ft.find(boundary_mapping["cylinder"])),
    V,
)

bcu = [bcu_inflow, bcu_walls, bcu_cylinder]

# Outlet
bcp_outlet = dirichletbc(
    PETSc.ScalarType(0),
    locate_dofs_topological(Q, fdim, ft.find(boundary_mapping["outlet"])),
    Q,
)
bcp = [bcp_outlet]

# --- Variational formulation --- #
u = TrialFunction(V)
v = TestFunction(V)
u_ = Function(V, name="u")
u_s = Function(V, name="u_tentative")
u_n = Function(V)
u_n1 = Function(V)
p = TrialFunction(Q)
q = TestFunction(Q)
p_ = Function(Q, name="p")
phi = Function(Q, name="phi")

# Define the first step, using Crank-Nicolson discretization, and semi-implicit
# Adams-Bashforth approximation. Uses three steps to solve for the fields.

f = Constant(mesh, PETSc.ScalarType((0, 0)))
F1 = rho / k * dot(u - u_n, v) * dx
F1 += inner(dot(1.5 * u_n - 0.5 * u_n1, 0.5 * nabla_grad(u + u_n)), v) * dx
F1 += 0.5 * mu * inner(grad(u + u_n), grad(v)) * dx - dot(p_, div(v)) * dx
F1 += dot(f, v) * dx
a1 = form(lhs(F1))
L1 = form(rhs(F1))
A1 = create_matrix(a1)
b1 = create_vector(extract_function_spaces(L1))

# Define the second step
a2 = form(dot(grad(p), grad(q)) * dx)
L2 = form(-rho / k * dot(div(u_s), q) * dx)
A2 = assemble_matrix(a2, bcs=bcp)
A2.assemble()
b2 = create_vector(extract_function_spaces(L2))

# Finally create the final step
a3 = form(rho * dot(u, v) * dx)
L3 = form(rho * dot(u_s, v) * dx - k * dot(nabla_grad(phi), v) * dx)
A3 = assemble_matrix(a3)
A3.assemble()
b3 = create_vector(extract_function_spaces(L3))


# Solver for step 1
solver1 = PETSc.KSP().create(mesh.comm)
solver1.setOperators(A1)
solver1.setType(PETSc.KSP.Type.BCGS)
pc1 = solver1.getPC()
pc1.setType(PETSc.PC.Type.JACOBI)

# Solver for step 2
solver2 = PETSc.KSP().create(mesh.comm)
solver2.setOperators(A2)
solver2.setType(PETSc.KSP.Type.MINRES)
pc2 = solver2.getPC()
pc2.setType(PETSc.PC.Type.HYPRE)
pc2.setHYPREType("boomeramg")

# Solver for step 3
solver3 = PETSc.KSP().create(mesh.comm)
solver3.setOperators(A3)
solver3.setType(PETSc.KSP.Type.CG)
pc3 = solver3.getPC()
pc3.setType(PETSc.PC.Type.SOR)


# --- Verification computing lift and drag coefficients --- #
n = -FacetNormal(mesh)
dObs = Measure(
    "ds", domain=mesh, subdomain_data=ft, subdomain_id=boundary_mapping["cylinder"]
)
u_t = inner(as_vector((n[1], -n[0])), u_)
# 2/(D*Ubar**2), not the tutorial's hardcoded `2/0.1` -- that constant is
# only equal to 0.1 when Ubar=1.0 (true for cases 2/3's Um=1.5, never
# generalized to case 1's Um=0.3 -> Ubar=0.2, D*Ubar**2=0.004).
norm_const = 2 / (D * Ubar**2)
drag = form(norm_const * (mu / rho * inner(grad(u_t), n) * n[1] - p_ * n[0]) * dObs)
lift = form(-norm_const * (mu / rho * inner(grad(u_t), n) * n[0] + p_ * n[1]) * dObs)
if mesh.comm.rank == 0:
    C_D = np.zeros(num_steps, dtype=PETSc.ScalarType)
    C_L = np.zeros(num_steps, dtype=PETSc.ScalarType)
    t_u = np.zeros(num_steps, dtype=np.float64)
    t_p = np.zeros(num_steps, dtype=np.float64)

tree = bb_tree(mesh, mesh.geometry.dim)
points = np.array([[0.15, 0.2, 0], [0.25, 0.2, 0]])
cell_candidates = compute_collisions_points(tree, points)
colliding_cells = compute_colliding_cells(mesh, cell_candidates, points)
front_cells = colliding_cells.links(0)
back_cells = colliding_cells.links(1)
if mesh.comm.rank == 0:
    p_diff = np.zeros(num_steps, dtype=PETSc.ScalarType)


folder = cylinderDir / "results" / "dokken" / meshFile.stem / str(args.inletBC)
folder.mkdir(exist_ok=True, parents=True)
vtx_u = VTXWriter(mesh.comm, folder / f"{case['name']}-u.bp", [u_], engine="BP4")
vtx_p = VTXWriter(mesh.comm, folder / f"{case['name']}-p.bp", [p_], engine="BP4")

# Write cadence: one frame every `write_every` steps, sized so a `fps`-fps
# animation over ALL `num_steps` frames plays back the whole `T`-second
# simulated interval (~num_steps/(fps*T) steps between frames) -- not tied to
# a fixed step count, so it stays correct if dt/T change. Previously wrote
# every step, which produced a 14G results/ dir for one 8s run (see
# cylinder/log.md, 2026-09-30 #fem entry).
animation_fps = 30
write_every = max(1, round(num_steps / (animation_fps * T)))
vtx_u.write(t)
vtx_p.write(t)
progress = tqdm.autonotebook.tqdm(desc="Solving PDE", total=num_steps)
for i in range(num_steps):
    progress.update(1)
    # Update current time step
    t += dt
    # Update inlet velocity
    inlet_velocity.t = t
    u_inlet.interpolate(inlet_velocity)

    # Step 1: Tentative velocity step
    A1.zeroEntries()
    assemble_matrix(A1, a1, bcs=bcu)
    A1.assemble()
    with b1.localForm() as loc:
        loc.set(0)
    assemble_vector(b1, L1)
    apply_lifting(b1, [a1], [bcu])
    b1.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    set_bc(b1, bcu)
    solver1.solve(b1, u_s.x.petsc_vec)
    u_s.x.scatter_forward()

    # Step 2: Pressure corrrection step
    with b2.localForm() as loc:
        loc.set(0)
    assemble_vector(b2, L2)
    apply_lifting(b2, [a2], [bcp])
    b2.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    set_bc(b2, bcp)
    solver2.solve(b2, phi.x.petsc_vec)
    phi.x.scatter_forward()

    p_.x.petsc_vec.axpy(1, phi.x.petsc_vec)
    p_.x.scatter_forward()

    # Step 3: Velocity correction step
    with b3.localForm() as loc:
        loc.set(0)
    assemble_vector(b3, L3)
    b3.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    solver3.solve(b3, u_.x.petsc_vec)
    u_.x.scatter_forward()

    # Write solutions to file -- every write_every steps only, not every
    # step (see write_every's definition above). Cd/Cl/p_diff below are
    # still recorded every step regardless (cheap scalars, not a bp write).
    if (i + 1) % write_every == 0 or i == num_steps - 1:
        vtx_u.write(t)
        vtx_p.write(t)

    # Update variable with solution form this time step
    with (
        u_.x.petsc_vec.localForm() as loc_,
        u_n.x.petsc_vec.localForm() as loc_n,
        u_n1.x.petsc_vec.localForm() as loc_n1,
    ):
        loc_n.copy(loc_n1)
        loc_.copy(loc_n)

    # Compute physical quantities
    # For this to work in paralell, we gather contributions from all processors
    # to processor zero and sum the contributions.
    drag_coeff = mesh.comm.gather(assemble_scalar(drag), root=0)
    lift_coeff = mesh.comm.gather(assemble_scalar(lift), root=0)
    p_front = None
    if len(front_cells) > 0:
        p_front = p_.eval(points[0], front_cells[:1])
    p_front = mesh.comm.gather(p_front, root=0)
    p_back = None
    if len(back_cells) > 0:
        p_back = p_.eval(points[1], back_cells[:1])
    p_back = mesh.comm.gather(p_back, root=0)
    if mesh.comm.rank == 0:
        t_u[i] = t
        t_p[i] = t - dt / 2
        C_D[i] = sum(drag_coeff)
        C_L[i] = sum(lift_coeff)
        # Choose first pressure that is found from the different processors
        for pressure in p_front:
            if pressure is not None:
                p_diff[i] = pressure[0]
                break
        for pressure in p_back:
            if pressure is not None:
                p_diff[i] -= pressure[0]
                break
progress.close()
vtx_u.close()
vtx_p.close()

A1.destroy()
A2.destroy()
A3.destroy()
b1.destroy()
b2.destroy()
b3.destroy()
solver1.destroy()
solver2.destroy()
solver3.destroy()

# Cd/Cl/p_diff were only ever held in memory and discarded at exit -- save
# them to their own csv (separate from cylinder/results.csv, which holds one
# summary row per case, not a per-timestep series). t_p != t_u (p_diff is
# evaluated half a step behind, per this scheme's pressure staggering), so
# both time columns are kept rather than assuming they line up.
if mesh.comm.rank == 0:
    coeffs_file = folder / "coefficients.csv"
    header = "t_u,C_D,C_L,t_p,p_diff"
    data = np.column_stack(
        [np.real(t_u), np.real(C_D), np.real(C_L), np.real(t_p), np.real(p_diff)]
    )
    np.savetxt(coeffs_file, data, delimiter=",", header=header, comments="")
    print(f"Saved Cd/Cl/p_diff to {coeffs_file}")

print("Script finished running!")
