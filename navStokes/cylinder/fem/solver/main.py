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
from solvers import Chorin
from yamlreader import read_config

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
bc_enf = config["solver"]["bc_enforcement"]
folder = results_dir / time_solver / bc_enf / caseName
folder.mkdir(exist_ok=True, parents=True)
writer = VTX_Files(folder)


# --- Velocity inlet boundary conditions ---
FlowClass = SinuPoiseuilleFlow if isRamped else ConstantPoiseuilleFlow
inlet_velocity = FlowClass(t=0.0, Um=Um)

chorin = Chorin(V, Q, faceTags, vDegree, writer)
chorin.set_physical_parameters(config)
u, p = chorin.solve_weak(dt, T, inlet_velocity)
writer.close()
