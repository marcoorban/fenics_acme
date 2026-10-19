import numpy as np 
import ufl 
from dolfinx import default_scalar_type, fem, geometry 
from dolfinx.fem.petsc import LinearProblem 
from dolfinx.io import XDMFFile, gmsh 
from mpi4py import MPI
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from pathlib import Path

# ---- READ MESH ---- #

fileName = "cavity-2D.msh"
meshFile = Path(fileName)
meshData = gmsh.read_from_msh(meshFile, MPI.COMM_WORLD, rank=0, gdim=2)

# ---- CREATE GEOMETRY ----# 
domain, cellTags, facetTags = ((meshData.mesh, meshData.cell_tags, meshData.facet_tags) 
                                if hasattr(meshData, "mesh")
                               else self.meshData)
tdim = domain.topology.dim 
fdim = domain.topology.dim -1 
domain.topology.create_connectivity(fdim, tdim)
n = ufl.FacetNormal(domain)
h = ufl.CellDiameter(domain)
dx = ufl.Measure("dx", domain=domain, subdomain_data=ct)
ds = ufl.Measure("ds", domain=domain, subdomain_data=ft)

# ---- CREATE FUNTION SPACE ----#
family = "Lagrange"
polyOrder = 2
V = fem.functionspace(domain, (family, polyOrder))
u = ufl.TrialFunction(V)
w = ufl.TestFunction(V)

# ----- WEAK FORM ------# 

# Compute the time derivative using first-order approximation 
du_dt = 
# First we code the bilinear form.
a = ufl.dot(w, ufl.p
# Now we represent the linear form.

# ----- SUPG ------- #


