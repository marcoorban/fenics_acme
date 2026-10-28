import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Column index (0 = coordinate) for each Reynolds number Ghia, Ghia &
# Shin (1982) actually tabulated -- see navStokes/cavity/benchmarks/.
GHIA_RE_COLUMNS = {100: 1, 400: 2, 1000: 3, 3200: 4, 5000: 5, 7500: 6, 10000: 7}

parser = argparse.ArgumentParser(
    description="Plot u along the vertical centerline and v along the "
    "horizontal centerline from a fem.py run, overlaid with the matching "
    "Ghia, Ghia & Shin (1982) benchmark columns."
)
parser.add_argument(
    "results_dir",
    help="Directory containing dat/u_vertical.dat and dat/v_horizontal.dat.",
)
parser.add_argument(
    "reynolds",
    type=float,
    help="Reynolds number this run solved for (for the Ghia overlay + lid-speed rescaling).",
)
args = parser.parse_args()

results_dir = Path(args.results_dir)
benchmarks_dir = Path(__file__).resolve().parent.parent / "benchmarks"

y, u = np.loadtxt(results_dir / "dat" / "u_vertical.dat", unpack=True)
x, v = np.loadtxt(results_dir / "dat" / "v_horizontal.dat", unpack=True)

nearest_ghia_re = min(GHIA_RE_COLUMNS, key=lambda re: abs(re - args.reynolds))
if nearest_ghia_re != args.reynolds:
    print(
        f"Note: Ghia et al. only tabulate Re={sorted(GHIA_RE_COLUMNS)}; "
        f"overlaying the closest one (Re={nearest_ghia_re}) against this "
        f"Re={args.reynolds:g} run."
    )
col = GHIA_RE_COLUMNS[nearest_ghia_re]
lid_speed = args.reynolds  # U_lid = Re*nu/length, with nu=length=1 in fem.py

ghia_y, ghia_u = np.loadtxt(
    benchmarks_dir / "ghia1982_table1_u_vertical.dat", usecols=(0, col), unpack=True
)
ghia_x, ghia_v = np.loadtxt(
    benchmarks_dir / "ghia1982_table2_v_horizontal.dat", usecols=(0, col), unpack=True
)

fig, axes = plt.subplots(1, 2, figsize=(11, 5))

axes[0].plot(u, y, color="#0072B2", lw=2, label="fem.py")
axes[0].plot(
    ghia_u * lid_speed,
    ghia_y,
    "o",
    color="black",
    ms=5,
    label=f"Ghia et al. (1982), Re={nearest_ghia_re}",
)
axes[0].set_xlabel("u")
axes[0].set_ylabel("y")
axes[0].set_title("Vertical centerline (x=0.5)")
axes[0].legend()
axes[0].grid(True, alpha=0.3)

axes[1].plot(x, v, color="#0072B2", lw=2, label="fem.py")
axes[1].plot(
    ghia_x,
    ghia_v * lid_speed,
    "o",
    color="black",
    ms=5,
    label=f"Ghia et al. (1982), Re={nearest_ghia_re}",
)
axes[1].set_xlabel("x")
axes[1].set_ylabel("v")
axes[1].set_title("Horizontal centerline (y=0.5)")
axes[1].legend()
axes[1].grid(True, alpha=0.3)

fig.suptitle(f"Re={args.reynolds:g}")
fig.tight_layout()

out_path = results_dir / "centerlines.png"
fig.savefig(out_path, dpi=150)
print(f"Wrote {out_path}")
