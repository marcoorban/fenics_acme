import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    description="Plot the pressure field and velocity streamlines from a "
    "fem.py run's saved field_grid.npz."
)
parser.add_argument(
    "results_dir",
    help="Directory containing field_grid.npz (e.g. results/Re_400/cavity_1968cells).",
)
args = parser.parse_args()

results_dir = Path(args.results_dir)
data = np.load(results_dir / "field_grid.npz")
x, y, u, v, p = data["x"], data["y"], data["u"], data["v"], data["p"]
re = float(data["Re"])

fig, ax = plt.subplots(figsize=(7, 6))

pressure = ax.contourf(x, y, p, levels=40, cmap="viridis")
fig.colorbar(pressure, ax=ax, label="pressure p")

# Default auto-seeded streamlines miss the corner vortices -- their
# speed is orders of magnitude below the bulk flow, so matplotlib's
# density-based seeding skips them. Draw the bulk flow first, then seed
# a few extra streamlines explicitly near all four corners.
ax.streamplot(x, y, u, v, color="0.15", density=1.2, linewidth=0.8)

corner_eps = 0.04
corner_seeds = np.array(
    [
        [corner_eps, corner_eps],
        [1 - corner_eps, corner_eps],
        [corner_eps, 1 - corner_eps],
        [1 - corner_eps, 1 - corner_eps],
    ]
)
ax.streamplot(
    x, y, u, v, start_points=corner_seeds, color="crimson", linewidth=0.8
)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(f"Lid-driven cavity, Re={re:g} -- pressure + velocity streamlines")

out_path = results_dir / "fields.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Wrote {out_path}")
