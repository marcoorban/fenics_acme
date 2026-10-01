This solver is based on the jsdokken.com dolfinx tutorials.
In particular, it uses as Crank-Nicolson discretization 
and a semi-implicit Adams-Bashroth approximation to 
numerical solve the Navier Stokes equation.

It also sets strong boundary conditions, and does not 
use any SUPG stabilization.
