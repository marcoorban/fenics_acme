import numpy as np
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
)


class Chorin:
    """Chorin's (1968) projection method -- the non-incremental fractional-step
    scheme: tentative velocity (no pressure), a pressure Poisson solve for the
    *full* new pressure (not a correction added to the old one, unlike Goda's
    method), then a velocity projection/correction step.

    Field naming matches fem/dokken/navStokes.py's IPCS (Goda's method)
    implementation where the role is the same, so the two solvers stay easy
    to compare side by side.
    """

    # Same {name: tag} convention as fem/dokken/navStokes.py and mesh/*.py --
    # shared by every mesh in cylinder/mesh/.
    BOUNDARY_MAPPING = {"inlet": 1, "outlet": 2, "walls": 3, "cylinder": 4}

    def __init__(self, V, Q, ft):
        self.V = V
        self.Q = Q
        self.mesh = V.mesh
        self.ft = ft  # facet tags -- needed to locate inlet/outlet/walls/cylinder

        self.u_n = Function(V)  # u^n -- current, known velocity
        self.u_s = Function(V, name="u_tentative")  # u* -- auxiliary/tentative velocity (step 1)
        self.u_ = Function(V, name="u")  # u^{n+1} -- corrected velocity (step 3)
        self.p_ = Function(Q, name="p")  # p^{n+1} -- solved directly each step, not incremented

    class _InletVelocity:
        """Same ramped parabolic profile as fem/dokken/navStokes.py's InletVelocity."""

        def __init__(self, t, Um):
            self.t = t
            self.Um = Um

        def set_time(self, t):
            """Advance to the inlet profile's velocity at time t."""
            self.t = t

        def __call__(self, x):
            values = np.zeros((2, x.shape[1]), dtype=PETSc.ScalarType)
            values[0] = 4 * self.Um * np.sin(self.t * np.pi / 8) * x[1] * (0.41 - x[1]) / (0.41**2)
            return values

    def solve_strong(self, dt, T, Re):
        """Solve from t=0 to T with strongly-enforced Dirichlet BCs (dokken's
        methodology: Crank-Nicolson diffusion, semi-implicit convection, same
        PETSc solver/preconditioner choices), using Chorin's splitting.
        Returns (u_, p_), the Function objects holding the final velocity and
        pressure fields.
        """
        V, Q, mesh, ft = self.V, self.Q, self.mesh, self.ft

        # Um/D/rho fixed at dokken's canonical DFG 2D-3 values; Re sweeps nu
        # (Ubar*D/Re) rather than Um -- see the caveat in the writeup above
        # this method if you wanted it the other way around.
        Um = 1.5
        D = 0.1
        rho_val = 1.0
        Ubar = (2 / 3) * Um
        nu_val = Ubar * D / Re

        t = 0.0
        num_steps = int(T / dt)
        k = Constant(mesh, PETSc.ScalarType(dt))
        rho = Constant(mesh, PETSc.ScalarType(rho_val))
        mu = Constant(mesh, PETSc.ScalarType(nu_val * rho_val))  # nu*rho, matches dokken's mu usage
        f = Constant(mesh, PETSc.ScalarType((0, 0)))

        fdim = mesh.topology.dim - 1

        # --- Strong boundary conditions (dokken's pattern) --- #
        u_inlet = Function(V)
        inlet_velocity = self._InletVelocity(t, Um)
        u_inlet.interpolate(inlet_velocity)
        bcu_inflow = dirichletbc(
            u_inlet, locate_dofs_topological(V, fdim, ft.find(self.BOUNDARY_MAPPING["inlet"]))
        )
        u_nonslip = np.array((0,) * mesh.geometry.dim, dtype=PETSc.ScalarType)
        bcu_walls = dirichletbc(
            u_nonslip, locate_dofs_topological(V, fdim, ft.find(self.BOUNDARY_MAPPING["walls"])), V
        )
        bcu_cylinder = dirichletbc(
            u_nonslip, locate_dofs_topological(V, fdim, ft.find(self.BOUNDARY_MAPPING["cylinder"])), V
        )
        bcu = [bcu_inflow, bcu_walls, bcu_cylinder]
        bcp_outlet = dirichletbc(
            PETSc.ScalarType(0),
            locate_dofs_topological(Q, fdim, ft.find(self.BOUNDARY_MAPPING["outlet"])),
            Q,
        )
        bcp = [bcp_outlet]

        # --- Step 1: tentative velocity -- dokken's F1, minus the pressure
        # term (that's the whole difference between Chorin and Goda here) --- #
        u = TrialFunction(V)
        v = TestFunction(V)
        F1 = rho / k * dot(u - self.u_n, v) * dx
        F1 += inner(dot(self.u_n, nabla_grad(0.5 * (u + self.u_n))), v) * dx
        F1 += 0.5 * mu * inner(grad(u + self.u_n), grad(v)) * dx
        # -f, not dokken's `+f` -- F1 is a "residual = 0" expression, and
        # rhs() only recovers the correctly-signed L(v)=<f,v> if f enters
        # the residual as -<f,v> (same sign logic as the pressure/viscous
        # terms above). Dormant in dokken's own tutorial since its f is
        # zero; confirmed wrong-signed otherwise against a Poisson check
        # (see cylinder/log.md).
        F1 -= dot(f, v) * dx
        a1 = form(lhs(F1))
        L1 = form(rhs(F1))
        A1 = create_matrix(a1)
        b1 = create_vector(extract_function_spaces(L1))

        # --- Step 2: pressure Poisson, solved directly for p_ (not a
        # correction phi added onto an old pressure, since step 1 above never
        # used a lagged pressure to begin with) --- #
        p = TrialFunction(Q)
        q = TestFunction(Q)
        a2 = form(dot(grad(p), grad(q)) * dx)
        L2 = form(-rho / k * dot(div(self.u_s), q) * dx)
        A2 = assemble_matrix(a2, bcs=bcp)
        A2.assemble()
        b2 = create_vector(extract_function_spaces(L2))

        # --- Step 3: velocity correction/projection, using p_'s gradient
        # directly (dokken's step 3 uses phi's gradient instead) --- #
        a3 = form(rho * dot(u, v) * dx)
        L3 = form(rho * dot(self.u_s, v) * dx - k * dot(nabla_grad(self.p_), v) * dx)
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
            b1.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
            set_bc(b1, bcu)
            solver1.solve(b1, self.u_s.x.petsc_vec)
            self.u_s.x.scatter_forward()

            # Step 2
            with b2.localForm() as loc:
                loc.set(0)
            assemble_vector(b2, L2)
            apply_lifting(b2, [a2], [bcp])
            b2.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
            set_bc(b2, bcp)
            solver2.solve(b2, self.p_.x.petsc_vec)
            self.p_.x.scatter_forward()

            # Step 3
            with b3.localForm() as loc:
                loc.set(0)
            assemble_vector(b3, L3)
            b3.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
            solver3.solve(b3, self.u_.x.petsc_vec)
            self.u_.x.scatter_forward()

            with self.u_.x.petsc_vec.localForm() as loc_, self.u_n.x.petsc_vec.localForm() as loc_n:
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
