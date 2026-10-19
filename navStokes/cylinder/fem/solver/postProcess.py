from pathlib import Path

import numpy as np
from dolfinx.fem import assemble_scalar, form
from dolfinx.geometry import bb_tree, compute_colliding_cells, compute_collisions_points
from dolfinx.io import VTXWriter
from ufl import FacetNormal, Measure, as_vector, grad, inner


class VTX_Files:
    def __init__(self, results_dir, time_solver, bc_enforcement, case_name, mesh_file):
        mesh_name = Path(mesh_file).stem
        self.folder = (
            Path(results_dir) / time_solver / bc_enforcement / case_name / mesh_name
        )
        self.folder.mkdir(parents=True, exist_ok=True)
        self.vtx_u = None
        self.vtx_p = None
        self._drag_form = None
        self._lift_form = None
        self._probe_points = None
        self._probe_cells = None
        self._coeffs_file = None
        self._cd_max = None
        self._cl_max = None
        self._cd_last = None
        self._cl_last = None
        self._dp_last = None

    def __call__(self, t, u_, p_):
        if self.vtx_u is None:
            comm = u_.function_space.mesh.comm
            self.vtx_u = VTXWriter(comm, self.folder / "u.bp", [u_], engine="BP4")
            self.vtx_p = VTXWriter(comm, self.folder / "p.bp", [p_], engine="BP4")
        self.vtx_u.write(t)
        self.vtx_p.write(t)

    def compute_coefficients(self, u_, p_, ft, cylinder_tag, mu, rho, D, Um):
        """Assemble and return (C_D, C_L) on the cylinder boundary (None on
        non-root ranks). Same DFG-benchmark formula already validated in
        fem/dokken/femSolver.py: normal flipped to point into the fluid,
        tangential velocity u_t = (n_y, -n_x).u_, norm_const = 2/(D*Ubar^2).
        Pure computation -- does not write anything; see write_coefficients.
        """
        mesh = u_.function_space.mesh
        if self._drag_form is None:
            n = -FacetNormal(mesh)
            dObs = Measure(
                "ds", domain=mesh, subdomain_data=ft, subdomain_id=cylinder_tag
            )
            u_t = inner(as_vector((n[1], -n[0])), u_)
            Ubar = (2 / 3) * Um
            norm_const = 2 / (D * Ubar**2)
            self._drag_form = form(
                norm_const * (mu / rho * inner(grad(u_t), n) * n[1] - p_ * n[0]) * dObs
            )
            self._lift_form = form(
                -norm_const
                * (mu / rho * inner(grad(u_t), n) * n[0] + p_ * n[1])
                * dObs
            )

        drag = mesh.comm.gather(assemble_scalar(self._drag_form), root=0)
        lift = mesh.comm.gather(assemble_scalar(self._lift_form), root=0)
        if mesh.comm.rank == 0:
            return float(np.real(sum(drag))), float(np.real(sum(lift)))
        return None

    def compute_pressure_diff(self, p_):
        """p(front) - p(back) at the two standard DFG-benchmark probe
        points -- (0.15, 0.2) just ahead of the cylinder, (0.25, 0.2) just
        behind it -- same points and gather-first-found-rank pattern as
        fem/dokken/femSolver.py. Collective (uses mesh.comm.gather): call
        this on every rank unconditionally, never behind a rank-dependent
        branch, or the gather deadlocks. Returns a float on rank 0 once
        both probes are located, else None (always None on non-root
        ranks, by mpi4py's gather convention).
        """
        mesh = p_.function_space.mesh
        if self._probe_cells is None:
            points = np.array([[0.15, 0.2, 0.0], [0.25, 0.2, 0.0]])
            tree = bb_tree(mesh, mesh.geometry.dim)
            cell_candidates = compute_collisions_points(tree, points)
            colliding_cells = compute_colliding_cells(mesh, cell_candidates, points)
            self._probe_points = points
            self._probe_cells = (colliding_cells.links(0), colliding_cells.links(1))

        front_cells, back_cells = self._probe_cells
        p_front = (
            p_.eval(self._probe_points[0], front_cells[:1])
            if len(front_cells) > 0
            else None
        )
        p_back = (
            p_.eval(self._probe_points[1], back_cells[:1])
            if len(back_cells) > 0
            else None
        )
        p_front = mesh.comm.gather(p_front, root=0)
        p_back = mesh.comm.gather(p_back, root=0)
        if mesh.comm.rank == 0:
            pf = next((v[0] for v in p_front if v is not None), None)
            pb = next((v[0] for v in p_back if v is not None), None)
            if pf is not None and pb is not None:
                return float(pf - pb)
        return None

    def write_coefficients(self, t, u_, p_, ft, cylinder_tag, mu, rho, D, Um):
        """Compute C_D/C_L/delta_p at time t and append a row to
        coefficients.dat, in the same folder as the VTX files. Writes
        immediately (not buffered to close()). Call this at the cadence
        you want recorded -- it carries no cadence logic of its own.
        """
        result = self.compute_coefficients(u_, p_, ft, cylinder_tag, mu, rho, D, Um)
        # Collective -- call on every rank before the rank-dependent branch
        # below, regardless of whether `result` ends up used.
        delta_p = self.compute_pressure_diff(p_)
        if result is None:
            return
        C_D, C_L = result
        self._cd_max = C_D if self._cd_max is None else max(self._cd_max, C_D)
        self._cl_max = C_L if self._cl_max is None else max(self._cl_max, C_L)
        self._cd_last, self._cl_last, self._dp_last = C_D, C_L, delta_p
        if self._coeffs_file is None:
            self._coeffs_file = open(self.folder / "coefficients.dat", "w")
            self._coeffs_file.write("t,C_D,C_L,delta_p\n")
        self._coeffs_file.write(f"{t},{C_D},{C_L},{delta_p}\n")
        self._coeffs_file.flush()

    def close(self):
        """Close the VTX/coefficients files and return this run's summary
        stats: (cd_max, cl_max, cd_last, cl_last, dp_last). Only valid on
        rank 0 (compute_coefficients/compute_pressure_diff only return a
        value there) and are None elsewhere.
        """
        if self.vtx_u is not None:
            self.vtx_u.close()
            self.vtx_p.close()
        if self._coeffs_file is not None:
            self._coeffs_file.close()
        return (
            self._cd_max,
            self._cl_max,
            self._cd_last,
            self._cl_last,
            self._dp_last,
        )
