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
kappa = 1 * 10**-6
kappa_xx = kappa
kappa_xy = 0
kappa_yx = 0
kappa_yy = kappa
K = ufl.as_tensor([[kappa_xx, kappa_xy], [kappa_yx, kappa_yy]])
alpha_deg = 63.4
alpha = np.deg2rad(alpha_deg)
beta_mag = 1.0
beta_x = beta_mag * np.cos(alpha)  # velocity x-component
beta_y = beta_mag * np.sin(alpha)  # velocity y-component
beta = fem.Constant(
    domain,
    np.array([beta_x, beta_y], dtype=default_scalar_type),
)

# Source term, coordinates, and boundary conditions
f = fem.Constant(domain, default_scalar_type(0.0))
x = ufl.SpatialCoordinate(domain)
# Boundary functions
#
strong_tags = [1, 2]
weak_tags = [3, 4, 5]

g_map = {
    3: fem.Constant(domain, default_scalar_type(0.0)),
    4: fem.Constant(domain, default_scalar_type(0.0)),
    5: fem.Constant(domain, default_scalar_type(0.0)),
}

# Check that g_map only contains weak tags
assert sorted(g_map) == weak_tags

# Strong boundary conditions
u_strong = fem.Constant(domain, default_scalar_type(1.0))
strong_facets = np.sort(np.concatenate([ft.find(t) for t in strong_tags]))
strong_dofs = fem.locate_dofs_topological(V, fdim, strong_facets)
bcs = [fem.dirichletbc(u_strong, strong_dofs, V)]

# Nitsche parameters
gamma = fem.Constant(domain, default_scalar_type(-1.0))
h = ufl.CellDiameter(domain)
n = ufl.FacetNormal(domain)
C = fem.Constant(domain, default_scalar_type(4.0))
knorm = ufl.dot(n, K * n)
ds_weak = ds(tuple(weak_tags))
inlet_tags = [5]
outlet_tags = [3, 4]
"""-------------- Bilinear terms -------------- """
a = -ufl.dot(ufl.grad(w), beta * u - K * ufl.grad(u)) * dx
a += w * (ufl.dot(-K * ufl.grad(u), n) + ufl.dot(beta, n * u)) * ds_weak
for tag in inlet_tags:
    a += (ufl.dot(-gamma * K * ufl.grad(w), n) - ufl.dot(beta, n * w)) * u * ds(tag)
for tag in outlet_tags:
    a += ufl.dot(-gamma * K * ufl.grad(w), n) * u * ds(tag)
a += C * knorm / h * w * u * ds_weak
"""--------------- Linear terms --------------- """
L = w * f * dx
for tag in inlet_tags:
    L += (
        (ufl.dot(-gamma * K * ufl.grad(w), n) - ufl.dot(beta, n * w))
        * g_map[tag]
        * ds(tag)
    )
for tag in outlet_tags:
    L += ufl.dot(-gamma * K * ufl.grad(w), n) * g_map[tag] * ds(tag)
for tag, g in g_map.items():
    L += C * knorm / h * w * g * ds(tag)

from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import XDMFFile

problem = LinearProblem(
    a,
    L,
    bcs=bcs,
    petsc_options_prefix="basic_linear_problem",
    petsc_options={
        "ksp_type": "preonly",
        "pc_type": "lu",
    },
)
u_h = problem.solve()
u_h.name = "u"

outname = "solutionWeakStrong"
with XDMFFile(domain.comm, f"{outname}.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(u_h)
