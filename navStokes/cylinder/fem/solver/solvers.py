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
    rhs,
    sqrt,
    sym,
)


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

    def set_physical_parameters(self, config, inlet_velocity):
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
        self.inlet_velocity = inlet_velocity
        self.Re = self.compute_Re(self.inlet_velocity.Um)

    def compute_Re(self, U):
        """Reynolds number for this solver's current nu/D, given the inlet
        peak velocity Um (Ubar = (2/3)*Um for this parabolic profile).
        Call after set_physical_parameters() so self.nuVal/self.D exist.
        """
        Ubar = (2 / 3) * U
        return Ubar * self.D / self.nuVal


class Chorin(NavStokesFEMSolver):
    """Chorin's (non-incremental) splitting scheme -- tentative velocity
    (Step 1), pressure Poisson (Step 2), velocity correction (Step 3) --
    with strongly-enforced (dokken-style) Dirichlet BCs on Step 1.

    ChorinWeak overrides boundary_conditions()/extra_step1_terms() below to
    switch Step 1 to Nitsche-enforced (weak) BCs with SUPG instead; the
    three-step split, Step 2/3, and the time-stepping loop are shared and
    live entirely in solve().
    """

    def __init__(self, V, Q, ft, vDegree, writer):
        super().__init__(V, Q, ft, vDegree, writer)
        self.u_n = Function(V, name="u_now")  # u^n; current, known velocity
        self.u_ = Function(V, name="u_next")  # u^(n+1); next velocity (next time step)
        self.u_s = Function(V, name="u_aux")  # Auxiliary / tentative velocity
        self.p_ = Function(Q, name="p_next")  # p^(n+1); next pressure field

    def boundary_conditions(self, dirichlet_boundaries, V, fdim):
        """Strong (dokken-style) Dirichlet enforcement: one dirichletbc per
        (tag, target) pair in dirichlet_boundaries. ChorinWeak overrides
        this to return [] instead -- BCs enforced weakly inside Step 1's
        F1 via extra_step1_terms() rather than on the matrix/space.
        """
        bcu = []
        for tag, target in dirichlet_boundaries:
            dofs = locate_dofs_topological(V, fdim, self.ft.find(tag))
            if isinstance(target, Function):
                bcu.append(dirichletbc(target, dofs))
            else:
                bcu.append(dirichletbc(target, dofs, V))
        return bcu

    def extra_step1_terms(self, dirichlet_boundaries, k, dt):
        """No extra terms for strong enforcement -- Step 1's F1 is just the
        base bulk terms built in solve(). ChorinWeak overrides this to add
        the Nitsche ibp/consistency/penalty terms plus SUPG.
        """
        return 0

    def solve(self, dt, T):
        """Solve from t=0 to T, dokken's Crank-Nicolson/semi-implicit
        methodology, Chorin's (non-incremental) splitting. Step 1's
        boundary-condition enforcement and any extra stabilization terms
        come from boundary_conditions()/extra_step1_terms() -- this class's
        versions give strong (dokken-style) Dirichlet enforcement with no
        extra terms; ChorinWeak overrides both for Nitsche + SUPG instead.
        Call set_physical_parameters(config) first -- this uses
        self.rhoVal/self.D/self.f from it, and overrides self.mu/self.nu to
        match inlet_velocity.Um (Re sweeps nu via Ubar*D/Re).
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

        t = 0.0
        num_steps = int(T / dt)
        k = Constant(mesh, PETSc.ScalarType(dt))
        fdim = mesh.topology.dim - 1

        # --- Inlet/walls/cylinder targets -- enforced as strong
        # dirichletbc()s or as Nitsche targets, depending on
        # boundary_conditions()'s override. --- #
        u_inlet = Function(V)
        u_inlet.interpolate(self.inlet_velocity)
        zero = Constant(mesh, PETSc.ScalarType((0, 0)))
        dirichlet_boundaries = [
            (self.BOUNDARY_MAPPING["inlet"], u_inlet),
            (self.BOUNDARY_MAPPING["walls"], zero),
            (self.BOUNDARY_MAPPING["cylinder"], zero),
        ]
        bcu = self.boundary_conditions(dirichlet_boundaries, V, fdim)

        bcp_outlet = dirichletbc(
            PETSc.ScalarType(0),
            locate_dofs_topological(Q, fdim, ft.find(self.BOUNDARY_MAPPING["outlet"])),
            Q,
        )
        bcp = [bcp_outlet]

        # --- Step 1: tentative velocity -- dokken's F1, minus the pressure
        # term (that's what makes this Chorin's method, not Goda's/IPCS),
        # plus whatever extra_step1_terms() adds (nothing here; Nitsche +
        # SUPG for ChorinWeak). --- #
        unsteady = self.rho / k * dot(u - self.u_n, w) * dx
        convective = inner(dot(self.u_n, nabla_grad(0.5 * (u + self.u_n))), w) * dx
        viscous = 0.5 * self.mu * inner(grad(u + self.u_n), grad(w)) * dx
        # -<f,w>, not dokken's `+<f,w>` -- F1 is a "residual = 0" expression,
        # and rhs() only recovers the correctly-signed L(w)=<f,w> if f enters
        # as -<f,w> (see cylinder/log.md, 2026-10-05 -- dokken's own tutorial
        # has this backwards, harmless there only because its f is zero).
        forcing = -dot(self.f, w) * dx
        F1 = unsteady + convective + viscous + forcing
        F1 = F1 + self.extra_step1_terms(dirichlet_boundaries, k, dt)
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

        # Steps between writes needed to hit self.fps output frames per
        # second of simulated time, given this run's dt.
        write_every = max(1, round(1 / (self.fps * dt)))

        progress = tqdm.autonotebook.tqdm(desc="Solving PDE", total=num_steps)
        for i in range(num_steps):
            progress.update(1)
            t += dt
            self.inlet_velocity.set_time(t)
            u_inlet.interpolate(self.inlet_velocity)

            # Step 1 -- bcu is [] for weak enforcement (the BCs are already
            # baked into F1/a1/L1 above via extra_step1_terms()), so these
            # calls are a no-op in that case and a real strong-BC
            # application otherwise.
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

            # Write solutions to file, and the lift/drag coefficients
            # alongside -- both gated by write_every (write_coefficients is
            # deliberately not called every step, only here).
            if (i + 1) % write_every == 0 or i == num_steps - 1:
                self.writer(t, self.u_, self.p_)
                self.writer.write_coefficients(
                    t,
                    self.u_,
                    self.p_,
                    self.ft,
                    self.BOUNDARY_MAPPING["cylinder"],
                    self.mu,
                    self.rho,
                    self.D,
                    self.inlet_velocity.Um,
                )

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


class ChorinWeak(Chorin):
    """Chorin with Nitsche-enforced (weak) velocity BCs and SUPG
    stabilization on Step 1, instead of strong dirichletbc()s. Step 2/3,
    the time loop, and the Re/nu/mu setup are all inherited unchanged from
    Chorin.solve() -- only the two hooks below differ.
    """

    def boundary_conditions(self, dirichlet_boundaries, V, fdim):
        return []

    def extra_step1_terms(self, dirichlet_boundaries, k, dt):
        """Nitsche (ibp/consistency/penalty) enforcement of
        dirichlet_boundaries, plus SUPG. The Nitsche terms use their own
        stress operator (2*nu*sym(grad(.))) -- a self-contained, standard
        Nitsche formulation independent of Step 1's CN-averaged bulk
        viscous term, not derived from it. SUPG's residual R1 *does* track
        the actual bulk terms from solve() (CN-averaged convection/
        viscous), since SUPG stabilizes whatever equation is really being
        solved.
        """
        u, w, mesh = self.u, self.w, self.mesh
        a = self.u_n  # frozen transport velocity for SUPG's tau/R1 below

        def sigma_u(uu):
            return 2 * self.nu * sym(grad(uu))

        n = FacetNormal(mesh)
        h = CellDiameter(mesh)
        beta = Constant(
            mesh, PETSc.ScalarType(10 * self.vDegree**2)
        )  # Nitsche penalty; tune if unstable

        ibp = 0
        consistency = 0
        penalty = 0
        for tag, g in dirichlet_boundaries:
            ds_tag = Measure(
                "ds", domain=mesh, subdomain_data=self.ft, subdomain_id=tag
            )
            ibp += -dot(w, dot(sigma_u(u), n)) * ds_tag
            consistency += -dot(dot(sigma_u(w), n), u - g) * ds_tag
            penalty += (beta / h) * dot(u - g, w) * ds_tag

        R1 = (
            self.rho / k * (u - self.u_n)
            + self.rho * dot(self.u_n, nabla_grad(0.5 * (u + self.u_n)))
            - 0.5 * self.mu * div(grad(u + self.u_n))
            - self.f
        )
        tau = 1 / sqrt(
            (2 / dt) ** 2
            + (2 * sqrt(dot(a, a)) / h) ** 2
            + (4 * self.nuVal / h**2) ** 2
        )
        SUPG = tau * dot(dot(a, nabla_grad(w)), R1) * dx

        return ibp + consistency + penalty + SUPG


class IPCS(NavStokesFEMSolver):
    pass
