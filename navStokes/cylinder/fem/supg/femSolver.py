import numpy as np

from argparser import parse_args
from meshreader import read_mesh
from yamlreader import read_config
from solvers import Chorin
from boundaryConditions import SinuPoiseuilleFlow, ConstantPoiseuilleFlow

from petsc4py import PETSc

from basix.ufl import element
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
    CellDiameter,
    FacetNormal,
    Identity,
    Measure,
    TestFunction,
    TrialFunction,
    as_vector,
    div,
    dot,
    dx,
    inner,
    outer,
    sym,
    lhs,
    grad,
    nabla_grad,
    rhs,
    sqrt,
)

from dolfinx.geometry import bb_tree, compute_collisions_points, compute_colliding_cells

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

# --- Physical parameters --- #

# --- Discretization parameters --- #
dt = config["discretization"]["dt"]
k = Constant(mesh, PETSc.ScalarType(dt))
function_spaces = config["discretization"]["function_spaces"]
vFamily = function_spaces["velocity"]["family"]
vDegree = function_spaces["velocity"]["degree"]
pFamily = function_spaces["pressure"]["family"]
pDegree = function_spaces["pressure"]["degree"]
num_steps = int(T/dt)

# --- Function spaces --- #
v_cg = element(vFamily, mesh.basix_cell(), vDegree, shape=(mesh.geometry.dim,))
s_cg = element(pFamily, mesh.basix_cell(), pDegree)
V = functionspace(mesh, v_cg) 
Q = functionspace(mesh, s_cg) 

# --- Mesh dimensions --- # 
fdim = mesh.topology.dim - 1
gdim = 2

chorin = Chorin(V, Q, faceTags)
u, p = chorin_weak.solve_weak(dt, T, )
