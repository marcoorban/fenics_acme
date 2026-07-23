"""
YAML configuration parsing for the advection-diffusion Nitsche solver.
"""

import yaml


def load_config(path: str = "config.yaml") -> dict:
    """Read the YAML parameter file and return as a dict."""
    with open(path, "r") as fh:
        cfg = yaml.safe_load(fh)
    return cfg
