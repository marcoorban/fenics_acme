import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import vtk
from matplotlib.collections import PolyCollection
from vtk.util.numpy_support import vtk_to_numpy

parser = argparse.ArgumentParser(
    description="Plot the local Peclet number Pe_h = |u|*h/(2*nu) over the "
    "domain, per mesh cell, from a saved velocity .vtu file -- shows where "
    "local convection dominates diffusion (Pe_h > 1), i.e. where a "
    "stabilization scheme like SUPG is actually doing work."
)
parser.add_argument("u_vtu", help="Path to the velocity .vtu file.")
parser.add_argument(
    "--nu",
    type=float,
    default=0.01,
    help="Kinematic viscosity (default: this project's nu=0.01).",
)
args = parser.parse_args()

reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(args.u_vtu)
reader.Update()
grid = reader.GetOutput()

# Cell-averaged velocity magnitude, from the point-wise "u" field -- same
# point-to-cell averaging postprocess.py uses for div(u).
point_to_cell = vtk.vtkPointDataToCellData()
point_to_cell.SetInputData(grid)
point_to_cell.Update()
u_cell = vtk_to_numpy(point_to_cell.GetOutput().GetCellData().GetArray("u"))
u_mag = np.linalg.norm(u_cell, axis=1)

points = vtk_to_numpy(grid.GetPoints().GetData())
n_cells = grid.GetNumberOfCells()
h = np.empty(n_cells)
corners = []
for i in range(n_cells):
    cell = grid.GetCell(i)
    ids = [cell.GetPointId(j) for j in range(cell.GetNumberOfPoints())]
    pts = points[ids]
    # Cell diameter: max pairwise distance among all of the cell's nodes
    # (corners + any higher-order edge nodes from the P2 velocity field,
    # which VTK stores on the same cells) -- a safe superset of the true
    # extent either way, and matches the spirit of UFL's CellDiameter
    # used in the solver's own SUPG term.
    diffs = pts[:, None, :2] - pts[None, :, :2]
    h[i] = np.sqrt((diffs**2).sum(axis=-1)).max()
    # Only the first 4 points are the cell's corners (VTK's quadratic-quad
    # point ordering: corners, then edge midpoints) -- the mesh geometry
    # is straight-sided even though the velocity field is P2, so the
    # corners alone draw the true cell outline for plotting.
    corners.append(pts[:4, :2])

nu = args.nu
pe_h = u_mag * h / (2 * nu)

fig, ax = plt.subplots(figsize=(7, 6))
poly = PolyCollection(corners, array=pe_h, cmap="viridis", edgecolors="none")
ax.add_collection(poly)
ax.set_xlim(points[:, 0].min(), points[:, 0].max())
ax.set_ylim(points[:, 1].min(), points[:, 1].max())
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(f"Local Peclet number  Pe_h = |u| h / (2 nu),  nu={nu:g}")
cbar = fig.colorbar(poly, ax=ax)
cbar.set_label("Pe_h")
fig.tight_layout()

out_path = Path(args.u_vtu).with_name("peclet_field.png")
fig.savefig(out_path, dpi=150)
print(f"Wrote Peclet field -> {out_path}")
print(f"Cells: {n_cells}")
print(f"Pe_h: min={pe_h.min():.3g} max={pe_h.max():.3g} mean={pe_h.mean():.3g}")
print(f"Fraction of cells with Pe_h > 1: {(pe_h > 1).mean() * 100:.1f}%")
