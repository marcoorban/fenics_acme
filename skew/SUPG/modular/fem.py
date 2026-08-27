import argparse
from pathlib import Path

import numpy as np
import project_io
import ufl
from dolfinx import default_scalar_type, fem, geometry
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import XDMFFile, gmsh
from mpi4py import MPI


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

    def create_physics(self):
        self.physics = Physics(self.params["physics"], self.domain)

    def create_supg(self):
        """Streamline operator and stabilisation parameter for SUPG.

        tau = (h / 2|beta|) * min(1, Pe / (3 p^2)), the piecewise-linear form
        of xi(Pe), with Pe built from the streamline diffusivity
        dot(beta, K beta) / dot(beta, beta) so that it stays correct if K
        becomes anisotropic. Built from the Constants rather than the raw
        floats so that changing physics.kappa.value at runtime is picked up
        without rebuilding.
        """
        p = self.physics
        order = self.polynomials["polyOrder"]
        beta_norm = ufl.sqrt(ufl.dot(p.beta, p.beta))
        kappa_beta = ufl.dot(p.beta, p.kappa * p.beta) / ufl.dot(p.beta, p.beta)
        Pe = beta_norm * self.h / (2.0 * kappa_beta)
        xi = ufl.min_value(1.0, Pe / (3.0 * order**2))
        self.tau = self.h / (2.0 * beta_norm) * xi
        self.supgL_w = ufl.dot(p.beta, ufl.grad(self.w))

    def allBCs(self):
        """Every boundary entry, however it ends up being imposed.

        StrongForm and WeakForm impose the whole boundary their own way and
        ignore the strong/weak split in the YAML; only StrongWeakForm reads
        the two keys separately.
        """
        bc = self.boundary_conditions
        return bc["strong"] + bc["weak"]

    def applyStrongBCs(self):
        """Boundary entries imposed by strong boundary conditions."""
        return []

    def applyWeakBCs(self):
        """Boundary entries imposed by adding Nitsche terms to the forms."""
        return

    def bilinearForm(self):
        u, w, p = self.u, self.w, self.physics
        flux = p.beta * u - p.kappa * ufl.grad(u)
        self.a = -ufl.dot(ufl.grad(w), flux) * self.dx
        self.a += self.supgL_w * self.tau * ufl.div(flux) * self.dx

    def linearForm(self):
        w, f, dx = self.w, self.physics.f, self.dx
        self.L = w * f * dx
        self.L += self.supgL_w * self.tau * f * dx

    def inflow(self):
        """Indicator that is 1 where the stream enters the domain, else 0.

        dot(beta, n) is negative where the stream runs against the outward
        normal, i.e. an inlet, and positive at an outlet. Evaluated per
        quadrature point, so a single tagged boundary may be partly inlet and
        partly outlet.
        """
        return ufl.conditional(
            ufl.lt(ufl.dot(self.physics.beta, self.n), 0.0), 1.0, 0.0
        )

    def outflow(self):
        """Indicator that is 1 where the stream leaves the domain, else 0."""
        return ufl.conditional(
            ufl.gt(ufl.dot(self.physics.beta, self.n), 0.0), 1.0, 0.0
        )

    def nitscheTerms(self, entry):
        """Bilinear and linear contributions imposing u = value on one tag.

        The volume form drops the boundary integral that integration by parts
        produces, so the diffusive flux -kappa*grad(u).n is restored here as
        the consistency term, together with the adjoint-consistency term
        (weighted by -gamma, the sign convention of nitsche_strong.py, so that
        gamma = -1 selects the non-symmetric variant) and the C*kappa/h
        penalty.

        The advective half of the flux is upwinded. Where the stream enters the
        domain the flux is carried by the prescribed value, so that term is
        masked by the inflow indicator and sits on the right-hand side. Where
        the stream leaves, the flux is carried by the unknown itself, so the
        outflow term restores the boundary integral dropped by the volume form
        and belongs on the left. Without it nothing constrains u on an outflow
        boundary once C*kappa/h is small, and the operator is near-singular
        there.
        """
        u, w, p = self.u, self.w, self.physics
        n, ds = self.n, self.ds(entry["tag"])
        gamma, C = self.nitsche["gamma"], self.nitsche["C"]
        g = fem.Constant(self.domain, default_scalar_type(entry["value"]))
        penalty = C * p.kappa / self.h
        beta_n = ufl.dot(p.beta, n)

        a = -p.kappa * ufl.dot(ufl.grad(u), n) * w * ds
        a += -gamma * p.kappa * ufl.dot(ufl.grad(w), n) * u * ds
        a += penalty * u * w * ds
        a += self.outflow() * beta_n * u * w * ds

        L = -gamma * p.kappa * ufl.dot(ufl.grad(w), n) * g * ds
        L += penalty * g * w * ds
        L += -self.inflow() * beta_n * g * w * ds

        return a, L

    def assemble(self):
        """Build the forms in the order the boundary treatment requires."""
        self.bilinearForm()
        self.linearForm()
        self.applyWeakBCs()

    def readMesh(self, meshFile):
        self.meshData = gmsh.read_from_msh(meshFile, MPI.COMM_WORLD, rank=0, gdim=2)

    def print_params(self):
        print(self.physics)
        print(self.boundary_conditions, self.nitsche, self.polynomials)

    def create_geometry(self):
        self.domain, self.ct, self.ft = (
            (self.meshData.mesh, self.meshData.cell_tags, self.meshData.facet_tags)
            if hasattr(self.meshData, "mesh")
            else self.meshData
        )
        self.tdim = self.domain.topology.dim
        self.fdim = self.tdim - 1
        self.domain.topology.create_connectivity(self.fdim, self.tdim)
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
        self.bcs = self.applyStrongBCs()
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

    def sample_line(self, y=0.5, n=201):
        """Evaluate u_h at n points along the line y=const.

        x spans the mesh's bounding box in x, independent of mesh
        resolution, so rows from different meshes/solvers line up
        point-for-point for a convergence comparison. Points that don't
        land in any local cell (relevant only for MPI runs, where each
        rank owns a partition of the mesh) are dropped.
        """
        coords = self.domain.geometry.x
        xmin, xmax = coords[:, 0].min(), coords[:, 0].max()
        xs = np.linspace(xmin, xmax, n)
        points = np.zeros((n, 3))
        points[:, 0] = xs
        points[:, 1] = y

        bb_tree = geometry.bb_tree(self.domain, self.domain.topology.dim)
        cell_candidates = geometry.compute_collisions_points(bb_tree, points)
        colliding_cells = geometry.compute_colliding_cells(self.domain, cell_candidates, points)

        found_points, cells = [], []
        for i in range(n):
            links = colliding_cells.links(i)
            if len(links) > 0:
                found_points.append(points[i])
                cells.append(links[0])

        found_points = np.array(found_points)
        u_values = self.u_h.eval(found_points, cells).flatten()
        return found_points[:, 0], u_values

    def write_results(self, outputFile):
        with XDMFFile(self.domain.comm, f"{outputFile}.xdmf", "w") as xdmf:
            xdmf.write_mesh(self.domain)
            xdmf.write_function(self.u_h)

    def write_line_sample(self, datFile, y=0.5, n=201):
        """Write x, u(x, y) samples to datFile, one file per run.

        The header records the cell count and solver class as metadata; the
        caller is expected to give each run (mesh resolution x boundary
        treatment) its own datFile, so a convergence study ends up as one
        file per case rather than one shared, ever-growing file.
        """
        x, u = self.sample_line(y=y, n=n)
        n_cells = self.domain.topology.index_map(self.tdim).size_global
        with open(datFile, "w") as f:
            f.write(f"# elements={n_cells} solver={type(self).__name__}\n")
            for xi, ui in zip(x, u):
                f.write(f"{xi:.16e} {ui:.16e}\n")

    def plot(self, plotFile):
        pass


class StrongForm(FEM_Solver):
    def applyStrongBCs(self):
        """Build the DirichletBC objects for the strongly imposed boundaries.

        Each entry gives a facet tag and a value; the tag is looked up in the
        facet MeshTags to get the facets, those become dofs of V, and the pair
        (value, dofs) becomes a DirichletBC. Returns the list for handing to
        LinearProblem.
        """
        bcs = []
        for bc in self.allBCs():
            facets = self.ft.find(bc["tag"])
            dofs = fem.locate_dofs_topological(self.V, self.fdim, facets)
            value = fem.Constant(self.domain, default_scalar_type(bc["value"]))
            bcs.append(fem.dirichletbc(value, dofs, self.V))
        return bcs


class WeakForm(FEM_Solver):
    def applyWeakBCs(self):
        """Impose every boundary weakly, adding its Nitsche terms to the forms."""
        for bc in self.allBCs():
            a, L = self.nitscheTerms(bc)
            self.a += a
            self.L += L


class StrongWeakForm(FEM_Solver):
    def applyStrongBCs(self):
        """Build the DirichletBC objects for the strongly imposed boundaries.

        Each entry gives a facet tag and a value; the tag is looked up in the
        facet MeshTags to get the facets, those become dofs of V, and the pair
        (value, dofs) becomes a DirichletBC. Returns the list for handing to
        LinearProblem.
        """
        bcs = []
        for bc in self.boundary_conditions["strong"]:
            facets = self.ft.find(bc["tag"])
            dofs = fem.locate_dofs_topological(self.V, self.fdim, facets)
            value = fem.Constant(self.domain, default_scalar_type(bc["value"]))
            bcs.append(fem.dirichletbc(value, dofs, self.V))
        return bcs

    def applyWeakBCs(self):
        """Impose only the boundaries tagged weak in the YAML via Nitsche."""
        for bc in self.boundary_conditions["weak"]:
            a, L = self.nitscheTerms(bc)
            self.a += a
            self.L += L


SOLVERS = {
    "strong": (StrongForm, "solutionStrong"),
    "weak": (WeakForm, "solutionWeak"),
    "weakStrong": (StrongWeakForm, "solutionStrongWeak"),
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Solve the skew advection-diffusion problem."
    )
    parser.add_argument("mesh_file", help="Path to the Gmsh .msh file")
    parser.add_argument("solver", choices=SOLVERS, help="Boundary condition treatment")
    parser.add_argument(
        "--sample-line",
        action="store_true",
        help=(
            "Append u(x, y=0.5) samples to output/line_samples.dat, for a "
            "convergence comparison across runs"
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    Solver, outputFile = SOLVERS[args.solver]
    outputDir = Path("output")
    outputDir.mkdir(exist_ok=True)
    outputFile = Path.joinpath(outputDir, f"{outputFile}_{Path(args.mesh_file).stem}")

    solver = Solver("params.yaml", args.mesh_file)
    solver.assemble()
    solver.solve()
    solver.write_results(outputFile)
    if args.sample_line:
        solver.write_line_sample(f"{outputFile}.dat")
    u = solver.u_h.x.array
    print(
        f"{Solver.__name__:16s} min={u.min():+.4f} max={u.max():+.4f} -> {outputFile}.xdmf"
    )
