"""
Advection-diffusion equation with weakly imposed Dirichlet BCs (Nitsche's method).
DOLFINx / FEniCSx API.

MONOLITHIC VERSION — everything in one file (config loading, mesh/physics
setup, and the Nitsche assembly/solve) for quick standalone testing of the
solver logic. The modularized version lives in advection_diffusion/
(config.py, solver.py, main.py, __init__.py).

All input parameters are read from config.yaml.
"""

import sys
import numpy as np
import yaml
from mpi4py import MPI
from petsc4py import PETSc

import ufl
from dolfinx import mesh, fem, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import XDMFFile


# ================================================================== #
#  Load configuration
# ================================================================== #
def load_config(path: str = "config.yaml") -> dict:
    """Read the YAML parameter file and return as a dict."""
    with open(path, "r") as fh:
        cfg = yaml.safe_load(fh)
    return cfg


# ================================================================== #
#  Build mesh from config
# ================================================================== #
def build_mesh(cfg: dict):
    """Create a rectangle mesh from the domain/mesh parameters."""
    dom = cfg["domain"]
    Lx, Ly = dom["Lx"], dom["Ly"]
    Nx, Ny = dom["Nx"], dom["Ny"]

    cell_map = {
        "triangle": mesh.CellType.triangle,
        "quadrilateral": mesh.CellType.quadrilateral,
    }
    cell_type = cell_map[dom["cell_type"]]

    domain = mesh.create_rectangle(
        MPI.COMM_WORLD,
        [np.array([0.0, 0.0]), np.array([Lx, Ly])],
        [Nx, Ny],
        cell_type,
    )
    return domain


# ================================================================== #
#  Build the function space
# ================================================================== #
def build_function_space(domain, cfg: dict):
    """Create the scalar H^1 function space from element parameters."""
    elem = cfg["element"]
    family = elem["family"]
    degree = elem["degree"]
    return fem.functionspace(domain, (family, degree))


# ================================================================== #
#  Build physical constants from config
# ================================================================== #
def build_physics(domain, cfg: dict):
    """
    Return the diffusivity tensor K (2x2 UFL matrix), the advection
    velocity beta as a DOLFINx Constant, and the diffusion/advection
    term switches (booleans) that turn each contribution on or off in
    the bilinear/linear forms.
    """
    phys = cfg["physics"]

    diffusion_on = phys.get("diffusion", True)
    advection_on = phys.get("advection", True)

    # --- diffusivity tensor ---
    d = phys["diffusivity"]
    K = ufl.as_tensor(
        [
            [default_scalar_type(d["kxx"]), default_scalar_type(d["kxy"])],
            [default_scalar_type(d["kyx"]), default_scalar_type(d["kyy"])],
        ]
    )

    # --- advection velocity ---
    vel = phys["velocity"]
    beta = fem.Constant(
        domain,
        np.array([vel["beta_x"], vel["beta_y"]], dtype=default_scalar_type),
    )

    return K, beta, diffusion_on, advection_on


# ================================================================== #
#  Mark boundary facets and create ds measure
# ================================================================== #
def build_boundary_measure(domain):
    """Tag all exterior facets with marker 1 and return ds measure."""
    tdim = domain.topology.dim
    fdim = tdim - 1
    domain.topology.create_connectivity(fdim, tdim)
    boundary_facets = mesh.exterior_facet_indices(domain.topology)

    facet_tags = mesh.meshtags(
        domain,
        fdim,
        boundary_facets,
        np.full(len(boundary_facets), 1, dtype=np.int32),
    )
    ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tags)
    inlet_facets = mesh.locate_entities_boundary(domain, fdim, inlet_marker)
    ds_inlet = ufl.Measure("ds", domain=domain, subdomain_data=inlet_facets)
    outlet_facets = mesh.locate_entities_boundary(domain, fdim, outlet_marker)
    ds_outlet = ufl.Measure("ds", domain=domain, subdomain_data=outlet_facets)
    return (ds, ds_inlet, ds_outlet)


def inlet_marker(x):
    """Locates the inlet face of a rectangular domain, assuming it is
    at x=0"""
    return np.isclose(x[0], 0.0)


def outlet_marker(x):
    """Located the outlet face of a rectangular domain, by reading
    the configuration file and obtaining the length of the rectangle."""
    return np.isclose(x[0], OUTLET)


# ================================================================== #
#  Assemble forms and solve
# ================================================================== #
def solve_nitsche(cfg: dict):
    """Top-level driver: build everything, solve, return solution."""

    # -- mesh & space ------------------------------------------------
    domain = build_mesh(cfg)
    V = build_function_space(domain, cfg)

    u = ufl.TrialFunction(V)
    w = ufl.TestFunction(V)

    # -- physics -----------------------------------------------------
    K, beta, diffusion_on, advection_on = build_physics(domain, cfg)

    # Source term (hardcoded for now — easy to extend)
    f = fem.Constant(domain, default_scalar_type(3.0))

    # Dirichlet data (hardcoded for now)
    x = ufl.SpatialCoordinate(domain)
    # g = ufl.sin(ufl.pi * x[1])
    g = x[0]

    # -- boundary measure --------------------------------------------
    (ds, ds_i, ds_o) = build_boundary_measure(domain)

    # -- Nitsche parameters ------------------------------------------
    # Check if gamma is either 1 or -1.
    gamma_coeff = cfg["nitsche"]["gamma"]
    if not (gamma_coeff == 1 or gamma_coeff == -1):
        sys.exit("Gamma must be 1 or -1!")
    gamma = fem.Constant(domain, default_scalar_type(cfg["nitsche"]["gamma"]))
    penalty = fem.Constant(domain, default_scalar_type(cfg["nitsche"]["penalty"]))
    h = ufl.CellDiameter(domain)
    n = ufl.FacetNormal(domain)

    # -- Bilinear form  a(u, v) -------------------------------------
    # a = 1 * u * v * ufl.dx  # zero form to start from (kept UFL-typed)

    # Volume: diffusion with full tensor K
    diffusion = ufl.inner(K * ufl.grad(u), ufl.grad(w)) * ufl.dx

    # Volume: advection
    advection = ufl.dot(beta, ufl.grad(u)) * w * ufl.dx

    if diffusion_on and advection_on:
        a = diffusion + advection
    elif diffusion_on and not advection_on:
        a = diffusion
    elif not diffusion_on and advection_on:
        a = advection
    else:
        a = None
        sys.exit(
            "Both diffusion and convection terms are turned off! Did you forget to turn them on?"
        )

    # -- Linear form  L(v) ------------------------------------------
    L = f * w * ufl.dx

    if diffusion_on:
        # Nitsche boundary (tag 1) — weakly imposed Dirichlet BC via the
        # diffusive flux. These terms only make sense when diffusion is on.
        # Normal flux through the tensor:  (K grad u) . n
        Kgradu_n = ufl.dot(-K * ufl.grad(u), n)
        Kgradv_n = ufl.dot(-gamma * K * ufl.grad(w), n)
        # Extra terms of bilinear form (LHS)
        a += Kgradu_n * w * ds(1)  # consistency
        a += Kgradv_n * u * ds(1)  # symmetry
        # Extra terms for linear form (RHS)
        L += Kgradv_n * g * ds(1)  # symmetry

        # Penalty terms
        C = 1
        breakpoint()
        knorm = ufl.sqrt(ufl.inner(K, K))
        # Penalty for bilinear form (LHS)
        a += (C * knorm / h) * w * u * ds(1)
        # Penalty for linear form (RHS)
        L += (C * knorm / h) * w * g * ds(1)

    if advection_on:
        a += w * ufl.dot(a, n) * u * ds(1)  # consistency
        a += ufl.dot((-1) * a, n) * w * u * ds_i(1)  # symmetry, bilinear
        L += ufl.dot((0) * a, n) * w * g * ds_i(1)  # symmetry, linear
    # -- Solve -------------------------------------------------------
    slv = cfg["solver"]
    problem = LinearProblem(
        a,
        L,
        bcs=[],
        petsc_options_prefix="basic_linear_problem",
        petsc_options={
            "ksp_type": slv["ksp_type"],
            "pc_type": slv["pc_type"],
        },
    )
    u_h = problem.solve()
    u_h.name = "u"

    # -- Diagnostics -------------------------------------------------
    error_form = fem.form((u_h - g) ** 2 * ds(1))
    error_local = fem.assemble_scalar(error_form)
    error_global = np.sqrt(domain.comm.allreduce(error_local, op=MPI.SUM))
    if domain.comm.rank == 0:
        print(f"L2 boundary error ||u_h - g||_dOmega = {error_global:.6e}")

    # -- Output ------------------------------------------------------
    outname = cfg["output"]["filename"]
    with XDMFFile(domain.comm, f"{outname}.xdmf", "w") as xdmf:
        xdmf.write_mesh(domain)
        xdmf.write_function(u_h)

    return u_h


# ================================================================== #
#  Entry point
# ================================================================== #
if __name__ == "__main__":
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config.yaml"
    cfg = load_config(config_path)
    OUTLET = cfg["domain"]["Lx"]
    solve_nitsche(cfg)
