import ufl
import project_io
import numpy as np
from mpi4py import MPI
from dolfinx import fem, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import XDMFFile, gmsh


class Physics:
    """Physical parameters of the advection-diffusion problem.

    Keeps the raw values read from the YAML alongside the dolfinx Constants
    the variational forms are written against, so a solver carries one
    attribute instead of one per coefficient.
    """

    def __init__(self, params, domain):
        self.kappa_value = params["kappa"]
        self.alpha = params["alpha"]
        self.beta_mag = params["beta_mag"]
        self.f_value = params["f"]
        self.create_constants(domain)

    def create_constants(self, domain):
        theta = np.deg2rad(self.alpha)
        direction = np.array([np.cos(theta), np.sin(theta)])
        self.kappa = fem.Constant(domain, default_scalar_type(self.kappa_value))
        self.beta = fem.Constant(domain, default_scalar_type(self.beta_mag * direction))
        self.f = fem.Constant(domain, default_scalar_type(self.f_value))

    def peclet(self, h):
        """Cell Peclet number for a cell size h."""
        return self.beta_mag * h / (2.0 * self.kappa_value)

    def __repr__(self):
        return (
            f"Physics(kappa={self.kappa_value}, alpha={self.alpha}, "
            f"beta_mag={self.beta_mag}, f={self.f_value})"
        )


class FEM_Solver:
    def __init__(self, paramsFile, meshFile):
        self.setParams(paramsFile)
        self.readMesh(meshFile)
        self.create_geometry()
        self.create_function_space()
        self.create_physics()
        self.create_supg()

    def setParams(self, paramsFile):
        params = project_io.readParams(paramsFile)
        self.params = params
        self.boundary_conditions = params["boundary_conditions"]
        self.nitsche = params["Nitsche"]
        self.polynomials = params["polynomials"]
        return

    def create_physics(self):
        self.physics = Physics(self.params["physics"], self.domain)
        return

    def create_supg(self):
        """Streamline operator and stabilisation parameter for SUPG.

        tau uses the doubly-asymptotic form tau = (h / 2|beta|) * xi(Pe),
        with xi(Pe) = coth(Pe) - 1/Pe and Pe the cell Peclet number. Built
        from the Constants rather than the raw floats so that changing
        physics.kappa.value at runtime is picked up without rebuilding.
        """
        p = self.physics
        beta_norm = ufl.sqrt(ufl.dot(p.beta, p.beta))
        Pe = beta_norm * self.h / (2.0 * p.kappa)
        xi = 1.0 / ufl.tanh(Pe) - 1.0 / Pe
        self.tau = self.h / (2.0 * beta_norm) * xi
        self.supgL_w = ufl.dot(p.beta, ufl.grad(self.w))
        return

    def readMesh(self, meshFile):
        self.meshData = gmsh.read_from_msh(meshFile, MPI.COMM_WORLD, rank=0, gdim=2)
        return

    def print_params(self):
        print(self.physics)
        print(self.boundary_conditions, self.nitsche, self.polynomials)

    def create_geometry(self):
        self.domain, self.ct, self.ft = (
            (self.meshData.mesh, self.meshData.cell_tags, self.meshData.facet_tags)
            if hasattr(self.meshData, "mesh")
            else self.meshData
        )
        tdim = self.domain.topology.dim
        fdim = tdim - 1
        self.domain.topology.create_connectivity(fdim, tdim)
        self.n = ufl.FacetNormal(self.domain)
        self.h = ufl.CellDiameter(self.domain)
        self.dx = ufl.Measure("dx", domain=self.domain, subdomain_data=self.ct)
        self.ds = ufl.Measure("ds", domain=self.domain, subdomain_data=self.ft)

    def create_function_space(self):
        family = self.polynomials["polyFamily"]
        order = self.polynomials["polyOrder"]
        self.V = fem.functionspace(self.domain, (family, order))
        self.u = ufl.TrialFunction(self.V)
        self.w = ufl.TestFunction(self.V)

    def solve(self):
        problem = LinearProblem(
            self.a,
            self.L,
            bcs=self.bcs,
            petsc_options_prefix="basic_linear_problem",
            petsc_options={"ksp_type": "preonly", "pc_type": "lu"},
        )
        self.u_h = problem.solve()
        self.u_h.name = "u"

    def postProcess(self):
        pass

    def write_results(self, outputFile):
        with XDMFFile(self.domain.comm, f"{outputFile}.xdmf", "w") as xdmf:
            xdmf.write_mesh(self.domain)
            xdmf.write_function(self.u_h)

    def plot(self, plotFile):
        pass


class StrongForm(FEM_Solver):
    def bilinearForm(self):
        u, w, p = self.u, self.w, self.physics
        self.a = -ufl.dot(ufl.grad(w), p.beta * u - p.kappa * ufl.grad(u)) * self.dx

    def linearForm(self):
        w, f, dx = self.w, self.physics.f, self.dx
        self.L = w * f * dx
        self.L += self.supgL_w * self.tau * f * dx


class WeakForm(FEM_Solver):
    def bilinearForm(self):
        pass

    def linearForm(self):
        pass


class StrongWeakForm(FEM_Solver):
    def bilinearForm(self):
        pass

    def linearForm(self):
        pass


if __name__ == "__main__":
    strong = StrongForm("params.yaml", "rectangle.msh")
    strong.print_params()
