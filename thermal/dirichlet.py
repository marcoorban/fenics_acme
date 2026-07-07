#!/usr/bin/env python
# coding: utf-8

# In[45]:


from pathlib import Path

from mpi4py import MPI
from petsc4py.PETSc import ScalarType  # type: ignore

import numpy as np

import ufl
from dolfinx import fem, io, mesh, plot
from dolfinx.fem.petsc import LinearProblem

import sys


# ## FEM setup and solution

# ### Problem parameters

# In[46]:


BL_x = 0.0
BL_y = 0.0
TR_x = 2.0
TR_y = 1.0
nx = int(sys.argv[1])
ny = int(sys.argv[2])
numElements = nx * ny
k = 1
A = 1
q = 1
beta = 0
P = 1
c = beta * P
T1 = 0
T2 = 1
Tinf = 0.5
LL = TR_x - BL_x

convergenceFile = Path(sys.argv[3])


def u_analytical(mod):
    return lambda x: (
        q * LL**2 / (2 * k) * (x[0] / LL - (x[0] / LL) ** 2)
        + (T2 - T1) * (x[0] / LL)
        + T1
    )


u_numpy = u_analytical(np)
u_ufl = u_analytical(ufl)


def create_mesh(nx, ny):
    msh = mesh.create_rectangle(
        comm=MPI.COMM_WORLD,
        points=((BL_x, BL_y), (TR_x, TR_y)),
        n=(nx, ny),
        cell_type=mesh.CellType.quadrilateral,
    )
    x = ufl.SpatialCoordinate(msh)
    return msh, x


def function_space(mesh, degree=1, ftype="Lagrange"):
    return fem.functionspace(mesh, (ftype, degree))


def dirichletbc(mesh, V):
    tdim = mesh.topology.dim
    fdim = tdim - 1
    bnd_left = mesh.locate_entities_boundary(
        mesh,
        dim=fdim,
        marker=lambda x: np.isclose(x[0], BL_x),
    )
    bnd_right = mesh.locate_entities_boundary(
        mesh,
        dim=fdim,
        marker=lambda x: np.isclose(x[0], TR_x),
    )
    dofs_left = fem.locate_dofs_topological(V=V, entity_dim=fdim, entities=bnd_left)
    dofs_right = fem.locate_dofs_topological(V=V, entity_dim=fdim, entities=bnd_right)

    # ### Dirichlet boundary conditions
    bc_left = fem.dirichletbc(value=ScalarType(T1), dofs=dofs_left, V=V)
    bc_right = fem.dirichletbc(value=ScalarType(T2), dofs=dofs_right, V=V)
    bcs = [bc_left, bc_right]
    return bcs


def forcing_function(msh):
    f0 = A * q + c * Tinf
    f = fem.Constant(msh, ScalarType(f0))
    return f


def neumannbc(msh):
    g = fem.Constant(msh, ScalarType(0))
    return g


def problem_setup(mesh, V):
    # Set up bilinear form
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    f = forcing_function(mesh)
    g = neumannbc(mesh)
    # Set up bilinear form
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    # Set up linear form
    L = ufl.inner(f, v) * ufl.dx + ufl.inner(g, v) * ufl.ds
    return a, L


def solve_system(a, L, bcs):
    problem = LinearProblem(
        a,
        L,
        bcs=bcs,
        petsc_options_prefix="demo_poisson_",
        petsc_options={
            "ksp_type": "preonly",
            "pc_type": "lu",
            "ksp_error_if_not_converged": True,
        },
    )

    return problem.solve()


def FEM():
    # read yaml
    mesh, x = create_mesh(nx, ny)
    V = function_space(mesh)
    # Dirichlet boundary conditions
    bcs_dirichlet = dirichletbc(mesh, V)
    # Set up problem
    a, L = problem_setup(mesh, V)
    # Solve and return solution
    fem_sln = solve_system(a, L, bcs_dirichlet)
    # Obtain analytical solution
    analytical = u_ufl(x)


def post_processing(fem_sln, analytical):
    secondary_vals()
    visualization()
    error_L2(fem_sln, analytical)
    error_inf(fem_sln, analytical)
    return


def output_results(msh, fem, analytical):
    out_folder = Path("dirichlet")
    out_folder.mkdir(parents=True, exist_ok=True)
    with io.XDMFFile(msh.comm, out_folder / f"FEMsln_{numElements}.xdmf", "w") as file:
        file.write_mesh(msh)
        file.write_function(fem)

    return


## Compute the L2 error
def error_L2(uh, u_ex, degree_raise=3):
    # Create higher order function space
    degree = uh.function_space.ufl_element().degree
    family = uh.function_space.ufl_element().family_name
    mesh = uh.function_space.mesh
    W = fem.functionspace(mesh, (family, degree + degree_raise))
    # Interpolate approximate solution
    u_W = fem.Function(W)
    u_W.interpolate(uh)

    # Interpolate exact solution, special handling if exact
    # solution is a ufl expression or a python lambda funtion
    u_ex_W = fem.Function(W)
    if isinstance(u_ex, ufl.core.expr.Expr):
        u_expr = fem.Expression(u_ex, W.element.interpolation_points)
        u_ex_W.interpolate(u_expr)
    else:
        u_ex_W.interpolate(u_ex)
    # Compute the error in the higher order function space
    e_W = fem.Function(W)
    e_W.x.array[:] = u_W.x.array - u_ex_W.x.array

    # Integrate the error
    error = fem.form(ufl.inner(e_W, e_W) * ufl.dx)
    error_local = fem.assemble_scalar(error)
    error_global = mesh.comm.allreduce(error_local, op=MPI.SUM)
    return np.sqrt(error_global)


def error_inf(uh, u_ex, degree_raise=3):
    return


def convergence(Ns, uh, u_ex):
    comm = uh.function_space.mesh.comm
    Es = error_L2(uh, u_ex)
    hs = 1.0 / np.sqrt(Ns)
    if comm.rank == 0:
        return hs, Es


l2_error = error_L2(uh, u_analytical)
hs, Es = convergence(numElements, uh, u_analytical)


with open(out_folder / convergenceFile, mode="a") as c:
    c.write(f"{hs}\t{Es}\n")
