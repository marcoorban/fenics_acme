import argparse
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import gmsh


def cell_diameters(mesh_file):
    """Return (centroids (M,2), diameters (M,)) for every 2D element in the
    mesh, matching UFL's CellDiameter: the max pairwise vertex distance
    within the cell (exact for a convex triangle or quad).

    Pure gmsh + numpy -- no dolfinx needed, so this runs fine in the
    lighter project venv instead of the fenics-env conda env.
    """
    gmsh.initialize()
    gmsh.open(mesh_file)
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    coords = {tag: node_coords[3 * i : 3 * i + 2] for i, tag in enumerate(node_tags)}

    centroids = []
    diameters = []
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2)
    for etype, tags, flat_nodes in zip(elem_types, elem_tags, elem_node_tags):
        _, _, _, num_nodes, *_ = gmsh.model.mesh.getElementProperties(etype)
        per_elem = np.asarray(flat_nodes).reshape(-1, num_nodes)
        for elem_nodes in per_elem:
            pts = np.array([coords[n] for n in elem_nodes])
            centroids.append(pts.mean(axis=0))
            diffs = pts[:, None, :] - pts[None, :, :]
            diameters.append(np.sqrt((diffs**2).sum(axis=-1)).max())
    gmsh.finalize()
    return np.array(centroids), np.array(diameters)


def nearest_value(points, centroids, values, chunk=2000):
    """Nearest-neighbor lookup of `values` (defined at `centroids`) onto
    `points`, chunked so the full pairwise-distance matrix never has to fit
    in memory at once (no scipy.spatial.cKDTree available here)."""
    out = np.empty(len(points))
    for i in range(0, len(points), chunk):
        block = points[i : i + chunk]
        d2 = ((block[:, None, :] - centroids[None, :, :]) ** 2).sum(axis=-1)
        out[i : i + chunk] = values[np.argmin(d2, axis=1)]
    return out


parser = argparse.ArgumentParser(
    description="Plot the local cell Peclet number Pe_K = |u|*h_K/(2*nu) "
    "over the domain, using the velocity field from a fem_*.py run's "
    "field_grid.npz and the actual mesh's per-cell diameter h_K. High Pe "
    "(rule of thumb: >~2) flags cells where convection dominates diffusion "
    "enough that stabilization (e.g. SUPG) could start to matter."
)
parser.add_argument("mesh", help="Path to the .msh file the run used.")
parser.add_argument(
    "results_dir", help="Directory containing that run's field_grid.npz."
)
parser.add_argument(
    "--nu", type=float, default=1.0, help="Kinematic viscosity (matches fem_*.py's nu=1.0 Constant)."
)
args = parser.parse_args()

results_dir = Path(args.results_dir)
data = np.load(results_dir / "field_grid.npz")
xs, ys, u, v, Re = data["x"], data["y"], data["u"], data["v"], float(data["Re"])

centroids, diameters = cell_diameters(args.mesh)

Xg, Yg = np.meshgrid(xs, ys)
grid_points = np.column_stack([Xg.ravel(), Yg.ravel()])
h_grid = nearest_value(grid_points, centroids, diameters).reshape(Xg.shape)

speed = np.sqrt(u**2 + v**2)
Pe = speed * h_grid / (2.0 * args.nu)

fig, ax = plt.subplots(figsize=(6, 5))
cf = ax.contourf(Xg, Yg, Pe, levels=40, cmap="inferno")
fig.colorbar(cf, ax=ax, label="Cell Peclet number")
# Rule-of-thumb instability threshold for unstabilized central-difference-like
# FEM convection -- cells above this are where SUPG would start to matter.
ax.contour(Xg, Yg, Pe, levels=[2.0], colors="cyan", linewidths=1.2, linestyles="--")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal")
ax.set_title(f"Cell Peclet number, Re={Re:g} (max={Pe.max():.2f}, mean={Pe.mean():.2f})")
fig.tight_layout()

out_path = results_dir / "peclet.png"
fig.savefig(out_path, dpi=150)
print(f"Wrote {out_path} (max Pe={Pe.max():.3f}, mean Pe={Pe.mean():.3f})")
