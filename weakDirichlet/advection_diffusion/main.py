"""
Command-line entry point for the advection-diffusion Nitsche solver.
"""

import sys

from .config import load_config
from .solver import solve_nitsche


def main(argv=None):
    """Load config (path optional, from argv[1]) and run the solver."""
    argv = sys.argv if argv is None else argv
    config_path = argv[1] if len(argv) > 1 else "config.yaml"

    cfg = load_config(config_path)
    return solve_nitsche(cfg)


if __name__ == "__main__":
    main()
