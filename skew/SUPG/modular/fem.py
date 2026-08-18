import ufl
import project_io
import numpy as np
from mpi4py import MPI
from dolfinx import fem, default_scalar_type
from dolfinx.io import gmsh


class FEM_Solver:
    def __init__(self, paramsFile, meshFile):
        self.setParams(paramsFile)
        self.readMesh(meshFile)
        self.create_geometry()
        self.create_function_space()

    def setParams(self, paramsFile):
        params = project_io.readParams(paramsFile)
        self.physics = params["physics"]
        self.boundary_conditions = params["boundary_conditions"]
        self.nitsche = params["Nitsche"]
        self.polynomials = params["polynomials"]
        return

    def readMesh(self, meshFile):
        self.meshData = gmsh.read_from_msh(meshFile, MPI.COMM_WORLD, rank=0, gdim=2)
        return

    def print_params(self):
        print(self.physics, self.boundary_conditions, self.nitsche, self.polynomials)

    def create_geometry(self):
        self.domain, self.ct, self.ft = (
            (self.meshData.mesh, self.meshData.cell_tags, self.meshData.facet_tags)
            if hasattr(self.meshData, "mesh")
            else self.meshData
        )
        tdim = self.domain.topology.dim
        fdim = tdim - 1
        self.domain.topology.create_connectivity(fdim, tdim)

    def create_function_space(self):
        family = self.polynomials["polyFamily"]
        order = self.polynomials["polyOrder"]
        V = fem.functionspace(self.domain, (family, order))
        self.u = ufl.TrialFunction(V)
        self.w = ufl.TestFunction(V)

    def solve(self):
        pass

    def postProcess(self):
        pass

    def write_results(self, outputFile):
        pass

    def plot(self, plotFile):
        pass


strong = FEM_Solver("params.yaml", "rectangle.msh")
strong.print_params()
