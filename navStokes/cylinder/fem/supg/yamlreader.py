from pathlib import Path

import yaml

DEFAULT_CONFIG = Path(__file__).parent / "config.yaml"


def read_config(path=DEFAULT_CONFIG):
    """Load config.yaml's physics/discretization/cases sections as a dict."""
    with open(path) as f:
        return yaml.safe_load(f)
