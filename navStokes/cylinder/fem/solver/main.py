import time

from basix.ufl import element
from dolfinx.fem import (
    Constant,
    functionspace,
)
from petsc4py import PETSc
from pathlib import Path
from argparser import parse_args
from boundaryConditions import ConstantPoiseuilleFlow, SinuPoiseuilleFlow
from meshreader import read_mesh
from postProcess import VTX_Files
from solvers import Chorin, ChorinWeak
from yamlreader import read_config

start_time = time.time()

args = parse_args()
config = read_config(args.yaml)
meshFile = args.mesh
(mesh, faceTags) = read_mesh(meshFile)

# --- Cases --- #
case = args.case
caseName = config["cases"][case]["name"]
Um = config["cases"][case]["Um"]
isRamped = config["cases"][case]["ramped"]
T = config["cases"][case]["T"]

# --- Discretization parameters --- #
dt = config["discretization"]["dt"]
k = Constant(mesh, PETSc.ScalarType(dt))
function_spaces = config["discretization"]["function_spaces"]
vFamily = function_spaces["velocity"]["family"]
vDegree = function_spaces["velocity"]["degree"]
pFamily = function_spaces["pressure"]["family"]
pDegree = function_spaces["pressure"]["degree"]
num_steps = int(T / dt)

# --- Function spaces --- #
v_cg = element(vFamily, mesh.basix_cell(), vDegree, shape=(mesh.geometry.dim,))
s_cg = element(pFamily, mesh.basix_cell(), pDegree)
V = functionspace(mesh, v_cg)
Q = functionspace(mesh, s_cg)

# --- Mesh dimensions --- #
fdim = mesh.topology.dim - 1
gdim = 2

# --- Prepare writer for saving pressure and velocity fields --
results_dir = Path(__file__).parent.parent.parent / "results"
time_solver = config["solver"]["time-solver"]
bc_enf = args.bc or config["solver"]["bc_enforcement"]
writer = VTX_Files(results_dir, time_solver, bc_enf, caseName, meshFile)


# --- Velocity inlet boundary conditions ---
FlowClass = SinuPoiseuilleFlow if isRamped else ConstantPoiseuilleFlow
inlet_velocity = FlowClass(t=0.0, Um=Um)

SolverClass = ChorinWeak if bc_enf == "weak" else Chorin
chorin = SolverClass(V, Q, faceTags, vDegree, writer)
chorin.set_physical_parameters(config, inlet_velocity)
u, p = chorin.solve(dt, T)
cd_max, cl_max, cd_last, cl_last, dp_last = writer.close()
wall_time = time.time() - start_time
if mesh.comm.rank == 0:
    num_cells = mesh.topology.index_map(mesh.topology.dim).size_global
    num_nodes = mesh.geometry.index_map().size_global
    print(
        f"{num_cells},{num_nodes},{wall_time},{cd_max},{cl_max},"
        f"{cd_last},{cl_last},{dp_last}"
    )
