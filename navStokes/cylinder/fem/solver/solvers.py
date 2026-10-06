import numpy as np
import tqdm.autonotebook
from dolfinx.fem import (
    Constant,
    Function,
    dirichletbc,
    extract_function_spaces,
    form,
    locate_dofs_topological,
)
from dolfinx.fem.petsc import (
    apply_lifting,
    assemble_matrix,
    assemble_vector,
    create_matrix,
    create_vector,
    set_bc,
)
from petsc4py import PETSc
from ufl import (
    CellDiameter,
    FacetNormal,
    Measure,
    TestFunction,
    TrialFunction,
    div,
    dot,
    dx,
    grad,
    inner,
    lhs,
    nabla_grad,
    outer,
    rhs,
    sqrt,
    sym,
)

from boundaryConditions import SinuPoiseuilleFlow


class NavStokesFEMSolver:
    # Same {name: tag} convention as fem/dokken/navStokes.py and mesh/*.py.
    BOUNDARY_MAPPING = {"inlet": 1, "outlet": 2, "walls": 3, "cylinder": 4}

    def __init__(self, V, Q, ft, vDegree, writer):
        self.V = V
        self.Q = Q
        self.u = TrialFunction(V)
        self.w = TestFunction(V)
        self.p = TrialFunction(Q)
        self.q = TestFunction(Q)
        self.mesh = V.mesh
        self.ft = ft
        self.vDegree = vDegree
        self.writer = writer

    def set_physical_parameters(self, config):
        self.muVal = config["physics"]["mu"]
        self.rhoVal = config["physics"]["rho"]
        self.nuVal = self.muVal / self.rhoVal
        self.mu = Constant(self.mesh, PETSc.ScalarType(self.muVal))
        self.rho = Constant(self.mesh, PETSc.ScalarType(self.rhoVal))
        self.nu = Constant(self.mesh, PETSc.ScalarType(self.nuVal))
        self.D = config["physics"]["D"]
        self.t = config["physics"]["t0"]
        self.fx = config["physics"]["fx"]
        self.fy = config["physics"]["fy"]
        self.f = Constant(self.mesh, PETSc.ScalarType((self.fx, self.fy)))
        self.fps = config["results"]["fps"]

    def compute_Re(self, Um):
        """Reynolds number for this solver's current nu/D, given the inlet
        peak velocity Um (Ubar = (2/3)*Um for this parabolic profile).
        Call after set_physical_parameters() so self.nuVal/self.D exist.
        """
        Ubar = (2 / 3) * Um
        return Ubar * self.D / self.nuVal


class Chorin(NavStokesFEMSolver):
    def __init__(self, V, Q, ft, vDegree, writer):
        super().__init__(V, Q, ft, vDegree, writer)
        self.u_n = Function(V, name="u_now")  # u^n; current, known velocity
        self.u_ = Function(V, name="u_next")  # u^(n+1); next velocity (next time step)
        self.u_s = Function(V, name="u_aux")  # Auxiliary / tentative velocity
        self.p_ = Function(Q, name="p_next")  # p^(n+1); next pressure field

    def solve_strong(self, dt, T, Re):
        """Solve from t=0 to T with strongly-enforced Dirichlet BCs, dokken's
        Crank-Nicolson/semi-implicit methodology, Chorin's (non-incremental)
        splitting. Call set_physical_parameters(config) first -- this uses
        self.rhoVal/self.D/self.f from it, and overrides self.mu/self.nu to
        match the Re passed in here (Um fixed at dokken's canonical DFG 2D-3
        value, 1.5; Re sweeps nu via Ubar*D/Re -- same assumption as
        fem/supg/chorin.py's solve_strong).
        Returns (u_, p_).
        """
        V, Q, mesh, ft, w, u, p, q = (
            self.V,
            self.Q,
            self.mesh,
            self.ft,
            self.w,
            self.u,
            self.p,
            self.q,
        )

        Um = 1.5
        Ubar = (2 / 3) * Um
        self.nuVal = Ubar * self.D / Re
        self.nu = Constant(mesh, PETSc.ScalarType(self.nuVal))
        self.mu = Constant(mesh, PETSc.ScalarType(self.nuVal * self.rhoVal))

        t = 0.0
        num_steps = int(T / dt)
        k = Constant(mesh, PETSc.ScalarType(dt))
        fdim = mesh.topology.dim - 1

        # --- Strong boundary conditions (dokken's pattern) --- #
        u_inlet = Function(V)
        inlet_velocity = SinuPoiseuilleFlow(t=t, Um=Um)
        u_inlet.interpolate(inlet_velocity)
        bcu_inflow = dirichletbc(
            u_inlet,
            locate_dofs_topological(V, fdim, ft.find(self.BOUNDARY_MAPPING["inlet"])),
        )
        u_nonslip = np.array((0,) * mesh.geometry.dim, dtype=PETSc.ScalarType)
        bcu_walls = dirichletbc(
            u_nonslip,
            locate_dofs_topological(V, fdim, ft.find(self.BOUNDARY_MAPPING["walls"])),
            V,
        )
        bcu_cylinder = dirichletbc(
            u_nonslip,
            locate_dofs_topological(
                V, fdim, ft.find(self.BOUNDARY_MAPPING["cylinder"])
            ),
            V,
        )
        bcu = [bcu_inflow, bcu_walls, bcu_cylinder]
        bcp_outlet = dirichletbc(
            PETSc.ScalarType(0),
            locate_dofs_topological(Q, fdim, ft.find(self.BOUNDARY_MAPPING["outlet"])),
            Q,
        )
        bcp = [bcp_outlet]

        # --- Step 1: tentative velocity -- dokken's F1, minus the pressure
        # term (that's what makes this Chorin's method, not Goda's/IPCS) --- #
        unsteady = self.rho / k * dot(u - self.u_n, w) * dx
        convective = inner(dot(self.u_n, nabla_grad(0.5 * (u + self.u_n))), w) * dx
        viscous = 0.5 * self.mu * inner(grad(u + self.u_n), grad(w)) * dx
        # -<f,w>, not dokken's `+<f,w>` -- F1 is a "residual = 0" expression,
        # and rhs() only recovers the correctly-signed L(w)=<f,w> if f enters
        # as -<f,w> (see cylinder/log.md, 2026-10-05 -- dokken's own tutorial
        # has this backwards, harmless there only because its f is zero).
        forcing = -dot(self.f, w) * dx
        F1 = unsteady + convective + viscous + forcing
        a1 = form(lhs(F1))
        L1 = form(rhs(F1))
        A1 = create_matrix(a1)
        b1 = create_vector(extract_function_spaces(L1))

        # --- Step 2: pressure Poisson, solved directly for p_ (Chorin: no
        # lagged pressure in step 1, so nothing to correct -- this solves
        # for the full new pressure, not an increment) --- #
        a2 = form(dot(grad(p), grad(q)) * dx)
        L2 = form(-self.rho / k * dot(div(self.u_s), q) * dx)
        A2 = assemble_matrix(a2, bcs=bcp)
        A2.assemble()
        b2 = create_vector(extract_function_spaces(L2))

        # --- Step 3: velocity correction/projection, using p_'s gradient
        # directly (Goda's step 3 would use a pressure-correction phi
        # instead) --- #
        a3 = form(self.rho * dot(u, w) * dx)
        L3 = form(
            self.rho * dot(self.u_s, w) * dx - k * dot(nabla_grad(self.p_), w) * dx
        )
        A3 = assemble_matrix(a3)
        A3.assemble()
        b3 = create_vector(extract_function_spaces(L3))

        solver1 = PETSc.KSP().create(mesh.comm)
        solver1.setOperators(A1)
        solver1.setType(PETSc.KSP.Type.BCGS)
        solver1.getPC().setType(PETSc.PC.Type.JACOBI)

        solver2 = PETSc.KSP().create(mesh.comm)
        solver2.setOperators(A2)
        solver2.setType(PETSc.KSP.Type.MINRES)
        pc2 = solver2.getPC()
        pc2.setType(PETSc.PC.Type.HYPRE)
        pc2.setHYPREType("boomeramg")

        solver3 = PETSc.KSP().create(mesh.comm)
        solver3.setOperators(A3)
        solver3.setType(PETSc.KSP.Type.CG)
        solver3.getPC().setType(PETSc.PC.Type.SOR)

        for i in range(num_steps):
            t += dt
            inlet_velocity.set_time(t)
            u_inlet.interpolate(inlet_velocity)

            # Step 1
            A1.zeroEntries()
            assemble_matrix(A1, a1, bcs=bcu)
            A1.assemble()
            with b1.localForm() as loc:
                loc.set(0)
            assemble_vector(b1, L1)
            apply_lifting(b1, [a1], [bcu])
            b1.ghostUpdate(
                addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE
            )
            set_bc(b1, bcu)
            solver1.solve(b1, self.u_s.x.petsc_vec)
            self.u_s.x.scatter_forward()

            # Step 2
            with b2.localForm() as loc:
                loc.set(0)
            assemble_vector(b2, L2)
            apply_lifting(b2, [a2], [bcp])
            b2.ghostUpdate(
                addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE
            )
            set_bc(b2, bcp)
            solver2.solve(b2, self.p_.x.petsc_vec)
            self.p_.x.scatter_forward()

            # Step 3
            with b3.localForm() as loc:
                loc.set(0)
            assemble_vector(b3, L3)
            b3.ghostUpdate(
                addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE
            )
            solver3.solve(b3, self.u_.x.petsc_vec)
            self.u_.x.scatter_forward()

            with (
                self.u_.x.petsc_vec.localForm() as loc_,
                self.u_n.x.petsc_vec.localForm() as loc_n,
            ):
                loc_.copy(loc_n)

        solver1.destroy()
        solver2.destroy()
        solver3.destroy()
        A1.destroy()
        A2.destroy()
        A3.destroy()
        b1.destroy()
        b2.destroy()
        b3.destroy()

        return self.u_, self.p_

    def solve_weak(self, dt, T, inlet_velocity):
        """Solve from t=0 to T with Nitsche-enforced (weak) velocity BCs and
        SUPG stabilization, using Chorin's splitting -- same 3 steps as
        solve_strong, except step 1's tentative-velocity system replaces the
        strong dirichletbc()s with Nitsche's method (ibp/consistency/penalty)
        and adds SUPG. Transport velocity for both convection and SUPG is
        self.u_n -- semi-implicit, same choice solve_strong already uses,
        keeps the system linear in u each step (settled 2026-10-05; the
        other options were a 2nd-order AB extrapolation like dokken's, which
        needs u_n1 back on Chorin, or a fully implicit a=u needing Newton).
        Step 1's stress operator has no pressure term (sigma_u, not the full
        sigma) -- Chorin's step 1 never has pressure, same restriction as
        solve_strong. Steps 2-3 are unchanged from solve_strong (outlet
        stays a strong/natural BC -- no reason to weakly enforce a condition
        that's already do-nothing). Call set_physical_parameters(config)
        first, same as solve_strong. Returns (u_, p_).
        """
        V, Q, mesh, ft, w, u, p, q = (
            self.V,
            self.Q,
            self.mesh,
            self.ft,
            self.w,
            self.u,
            self.p,
            self.q,
        )

        Um = inlet_velocity.Um
        Ubar = (2 / 3) * Um
        Re = self.compute_Re(Um)
        self.nuVal = Ubar * self.D / Re
        self.nu = Constant(mesh, PETSc.ScalarType(self.nuVal))
        self.mu = Constant(mesh, PETSc.ScalarType(self.nuVal * self.rhoVal))

        t = 0.0
        num_steps = int(T / dt)
        k = Constant(mesh, PETSc.ScalarType(dt))
        fdim = mesh.topology.dim - 1

        a = self.u_n  # semi-implicit transport velocity

        # --- Inlet velocity -- still time-dependent, but now its value is a
        # Nitsche target `g` (consistency/penalty terms below) instead of a
        # dirichletbc. --- #
        u_inlet = Function(V)
        u_inlet.interpolate(inlet_velocity)
        zero = Constant(mesh, PETSc.ScalarType((0, 0)))
        dirichlet_boundaries = [
            (self.BOUNDARY_MAPPING["inlet"], u_inlet),
            (self.BOUNDARY_MAPPING["walls"], zero),
            (self.BOUNDARY_MAPPING["cylinder"], zero),
        ]

        # Outlet pressure BC stays strong, same as solve_strong.
        bcp_outlet = dirichletbc(
            PETSc.ScalarType(0),
            locate_dofs_topological(Q, fdim, ft.find(self.BOUNDARY_MAPPING["outlet"])),
            Q,
        )
        bcp = [bcp_outlet]

        # --- Step 1: tentative velocity, Nitsche BCs + SUPG --- #
        def sigma_u(uu):
            return 2 * self.nu * sym(grad(uu))

        n = FacetNormal(mesh)
        h = CellDiameter(mesh)
        beta = Constant(
            mesh, PETSc.ScalarType(10 * self.vDegree**2)
        )  # Nitsche penalty; tune if unstable

        unsteady = self.rho / k * dot(u - self.u_n, w) * dx
        convection = self.rho * inner(nabla_grad(w), outer(a, u)) * dx
        viscous = inner(sigma_u(u), grad(w)) * dx
        forcing = -dot(w, self.f) * dx

        ibp = 0
        consistency = 0
        penalty = 0
        for tag, g in dirichlet_boundaries:
            ds_tag = Measure("ds", domain=mesh, subdomain_data=ft, subdomain_id=tag)
            ibp += -dot(w, dot(sigma_u(u), n)) * ds_tag
            consistency += -dot(dot(sigma_u(w), n), u - g) * ds_tag
            penalty += (beta / h) * dot(u - g, w) * ds_tag

        R1 = (
            self.rho / k * (u - self.u_n)
            + self.rho * dot(a, nabla_grad(u))
            - div(sigma_u(u))
            - self.f
        )
        tau = 1 / sqrt(
            (2 / dt) ** 2
            + (2 * sqrt(dot(a, a)) / h) ** 2
            + (4 * self.nuVal / h**2) ** 2
        )
        SUPG = tau * dot(dot(a, nabla_grad(w)), R1) * dx

        F1 = (
            unsteady
            + convection
            + viscous
            + forcing
            + ibp
            + consistency
            + penalty
            + SUPG
        )
        a1 = form(lhs(F1))
        L1 = form(rhs(F1))
        A1 = create_matrix(a1)
        b1 = create_vector(extract_function_spaces(L1))

        # --- Step 2: pressure Poisson -- unchanged from solve_strong --- #
        a2 = form(dot(grad(p), grad(q)) * dx)
        L2 = form(-self.rho / k * dot(div(self.u_s), q) * dx)
        A2 = assemble_matrix(a2, bcs=bcp)
        A2.assemble()
        b2 = create_vector(extract_function_spaces(L2))

        # --- Step 3: velocity correction -- unchanged from solve_strong --- #
        a3 = form(self.rho * dot(u, w) * dx)
        L3 = form(
            self.rho * dot(self.u_s, w) * dx - k * dot(nabla_grad(self.p_), w) * dx
        )
        A3 = assemble_matrix(a3)
        A3.assemble()
        b3 = create_vector(extract_function_spaces(L3))

        solver1 = PETSc.KSP().create(mesh.comm)
        solver1.setOperators(A1)
        solver1.setType(PETSc.KSP.Type.BCGS)
        solver1.getPC().setType(PETSc.PC.Type.JACOBI)

        solver2 = PETSc.KSP().create(mesh.comm)
        solver2.setOperators(A2)
        solver2.setType(PETSc.KSP.Type.MINRES)
        pc2 = solver2.getPC()
        pc2.setType(PETSc.PC.Type.HYPRE)
        pc2.setHYPREType("boomeramg")

        solver3 = PETSc.KSP().create(mesh.comm)
        solver3.setOperators(A3)
        solver3.setType(PETSc.KSP.Type.CG)
        solver3.getPC().setType(PETSc.PC.Type.SOR)

        # Steps between writes needed to hit self.fps output frames per
        # second of simulated time, given this run's dt.
        write_every = max(1, round(1 / (self.fps * dt)))

        progress = tqdm.autonotebook.tqdm(desc="Solving PDE", total=num_steps)
        for i in range(num_steps):
            progress.update(1)
            t += dt
            inlet_velocity.set_time(t)
            u_inlet.interpolate(inlet_velocity)

            # Step 1 -- no bcs=[...] here, unlike solve_strong: the BCs are
            # inside F1 itself (Nitsche), not enforced on the matrix/space.
            A1.zeroEntries()
            assemble_matrix(A1, a1)
            A1.assemble()
            with b1.localForm() as loc:
                loc.set(0)
            assemble_vector(b1, L1)
            b1.ghostUpdate(
                addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE
            )
            solver1.solve(b1, self.u_s.x.petsc_vec)
            self.u_s.x.scatter_forward()

            # Step 2
            with b2.localForm() as loc:
                loc.set(0)
            assemble_vector(b2, L2)
            apply_lifting(b2, [a2], [bcp])
            b2.ghostUpdate(
                addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE
            )
            set_bc(b2, bcp)
            solver2.solve(b2, self.p_.x.petsc_vec)
            self.p_.x.scatter_forward()

            # Step 3
            with b3.localForm() as loc:
                loc.set(0)
            assemble_vector(b3, L3)
            b3.ghostUpdate(
                addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE
            )
            solver3.solve(b3, self.u_.x.petsc_vec)
            self.u_.x.scatter_forward()

            with (
                self.u_.x.petsc_vec.localForm() as loc_,
                self.u_n.x.petsc_vec.localForm() as loc_n,
            ):
                loc_.copy(loc_n)

            # Write solutions to file
            if (i + 1) % write_every == 0 or i == num_steps - 1:
                self.writer(t, self.u_, self.p_)

        progress.close()
        solver1.destroy()
        solver2.destroy()
        solver3.destroy()
        A1.destroy()
        A2.destroy()
        A3.destroy()
        b1.destroy()
        b2.destroy()
        b3.destroy()

        return self.u_, self.p_
