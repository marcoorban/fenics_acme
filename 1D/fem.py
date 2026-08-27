import argparse
from pathlib import Path

import numpy as np
import project_io
import ufl
from dolfinx import default_scalar_type, fem, geometry
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import XDMFFile, gmsh
from mpi4py import MPI
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
