from dolfinx.io import gmsh
from mpi4py import MPI

out = gmsh.read_from_msh("rectangle.msh", MPI.COMM_WORLD, rank=0, gdim=2)
domain, ct, ft = (
    (out.mesh, out.cell_tags, out.facet_tags) if hasattr(out, "mesh") else out
)

import ufl

tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)

ds = ufl.Measure("ds", domain=domain, subdomain_data=ft)
dx = ufl.Measure("dx", domain=domain, subdomain_data=ct)

from dolfinx import fem, default_scalar_type
import numpy as np
import ufl

V = fem.functionspace(domain, ("Lagrange", 1))
u = ufl.TrialFunction(V)
w = ufl.TestFunction(V)

# Physics parameters
kappa_xx = 0.01
kappa_xy = 0
kappa_yx = 0
kappa_yy = 0.01
K = ufl.as_tensor([[kappa_xx, kappa_xy], [kappa_yx, kappa_yy]])
beta_x = 1.0  # velocity x-component
beta_y = 0.0  # velocity y-component
beta = fem.Constant(
    domain,
    np.array([beta_x, beta_y], dtype=default_scalar_type),
)

# Source term, coordinates, and boundary conditions
f = fem.Constant(domain, default_scalar_type(0.0))
x = ufl.SpatialCoordinate(domain)
# Boundary functions
g_map = {
    1: fem.Constant(domain, default_scalar_type(1.0)),
    3: fem.Constant(domain, default_scalar_type(0.0)),
}

# Nitsche parameters
gamma = fem.Constant(domain, default_scalar_type(-1.0))
h = ufl.CellDiameter(domain)
n = ufl.FacetNormal(domain)
C = fem.Constant(domain, default_scalar_type(4.0))
knorm = ufl.dot(n, K * n)
inlet_tag = 1
outlet_tag = 3
"""-------------- Bilinear terms -------------- """
a = -ufl.dot(ufl.grad(w), beta * u - K * ufl.grad(u)) * dx
a += w * (ufl.dot(-K * ufl.grad(u), n) + ufl.dot(beta, n * u)) * ds
a += (ufl.dot(-gamma * K * ufl.grad(w), n) - ufl.dot(beta, n * w)) * u * ds(inlet_tag)
a += ufl.dot(-gamma * K * ufl.grad(w), n) * u * ds(outlet_tag)
a += C * knorm / h * w * u * ds
"""--------------- Linear terms --------------- """
L = w * f * dx
L += (
    (ufl.dot(-gamma * K * ufl.grad(w), n) - ufl.dot(beta, n * w))
    * g_map[inlet_tag]
    * ds(inlet_tag)
)
L += ufl.dot(-gamma * K * ufl.grad(w), n) * g_map[outlet_tag] * ds(outlet_tag)
for tag, g in g_map.items():
    L += C * knorm / h * w * g * ds(tag)

from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import XDMFFile

problem = LinearProblem(
    a,
    L,
    bcs=[],
    petsc_options_prefix="basic_linear_problem",
    petsc_options={
        "ksp_type": "preonly",
        "pc_type": "lu",
    },
)
u_h = problem.solve()
u_h.name = "u"

outname = "solution"
with XDMFFile(domain.comm, f"{outname}.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(u_h)
