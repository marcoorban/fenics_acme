"""
advection_diffusion
====================

Package for solving the advection-diffusion equation with weakly imposed
Dirichlet boundary conditions (Nitsche's method) using DOLFINx / FEniCSx.

Modules
-------
config
    YAML configuration loading.
solver
    Mesh/function-space construction and the Nitsche assembly/solve routine.
main
    Command-line entry point.
"""

from .config import load_config
from .solver import (
    build_mesh,
    build_function_space,
    build_physics,
    build_boundary_measure,
    solve_nitsche,
)
from .main import main

__all__ = [
    "load_config",
    "build_mesh",
    "build_function_space",
    "build_physics",
    "build_boundary_measure",
    "solve_nitsche",
    "main",
]
