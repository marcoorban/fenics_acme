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
print(meshData)
