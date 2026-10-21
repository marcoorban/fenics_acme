import argparse
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

parser = argparse.ArgumentParser(
    description="Sample u, v and div(u) along the cavity's horizontal "
    "(y=0.5) AND vertical (x=0.5) centerlines at 129 points each, for the "
    "strong/weak/weak_noSUPG velocity .vtu files of one mesh+Re case, and "
    "write them into two gnuplot-friendly .dat files under results/dat/ "
    "-- one per line, since Ghia et al. (1982) tabulates u along the "
    "vertical line and v along the horizontal one, not both along either."
)
parser.add_argument("strong_vtu", help="Path to fem_strong.py's velocity .vtu file.")
parser.add_argument("weak_vtu", help="Path to fem_weak.py's velocity .vtu file.")
parser.add_argument("weak_nosupg_vtu", help="Path to fem_weak_noSUPG.py's velocity .vtu file.")
args = parser.parse_args()

n_points = 129


def make_line(point1, point2):
    line = vtk.vtkLineSource()
    line.SetPoint1(*point1)
    line.SetPoint2(*point2)
    line.SetResolution(n_points - 1)
    line.Update()
    return line


def sample(vtu_path, line, coord_index):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_path)
    reader.Update()
    grid = reader.GetOutput()

    # Point-wise divergence of the velocity field, added as a point array
    # alongside "u" -- same gradient computation as postprocess.py, but
    # left per-point (not averaged to cells) so it can be probed at the
    # line's sample points like any other point field.
    gradient = vtk.vtkGradientFilter()
    gradient.SetInputData(grid)
    gradient.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "u")
    gradient.SetComputeDivergence(True)
    gradient.SetDivergenceArrayName("div_u")
    gradient.Update()

    probe = vtk.vtkProbeFilter()
    probe.SetInputConnection(line.GetOutputPort())
    probe.SetSourceData(gradient.GetOutput())
    probe.Update()
    sampled = probe.GetOutput()

    coord = vtk_to_numpy(sampled.GetPoints().GetData())[:, coord_index]
    velocity = vtk_to_numpy(sampled.GetPointData().GetArray("u"))
    div_u = vtk_to_numpy(sampled.GetPointData().GetArray("div_u"))
    return coord, velocity[:, 0], velocity[:, 1], div_u


def write_dat(out_path, coord_name, coord, strong, weak, nosupg):
    u_strong, v_strong, div_strong = strong
    u_weak, v_weak, div_weak = weak
    u_nosupg, v_nosupg, div_nosupg = nosupg
    np.savetxt(
        out_path,
        np.column_stack(
            [
                coord,
                u_strong, v_strong, div_strong,
                u_weak, v_weak, div_weak,
                u_nosupg, v_nosupg, div_nosupg,
            ]
        ),
        header=f"{coord_name} u_strong v_strong div_strong u_weak v_weak "
        "div_weak u_noSUPG v_noSUPG div_noSUPG",
        comments="# ",
    )
    print(f"Wrote centerline profile -> {out_path}")


# All three runs share the same mesh and Re, so both are read off the
# strong-BC path (results/strong/<mesh_stem>/Re_<re>/u_p0_000000.vtu) to
# name the output files.
mesh_stem = Path(args.strong_vtu).parts[-3]
re_label = Path(args.strong_vtu).parts[-2]
out_dir = Path(__file__).parent / "results" / "dat"
out_dir.mkdir(parents=True, exist_ok=True)

horizontal = make_line((0.0, 0.5, 0.0), (1.0, 0.5, 0.0))
x, *strong_h = sample(args.strong_vtu, horizontal, 0)
_, *weak_h = sample(args.weak_vtu, horizontal, 0)
_, *nosupg_h = sample(args.weak_nosupg_vtu, horizontal, 0)
write_dat(out_dir / f"{mesh_stem}_{re_label}.dat", "x", x, strong_h, weak_h, nosupg_h)

vertical = make_line((0.5, 0.0, 0.0), (0.5, 1.0, 0.0))
y, *strong_v = sample(args.strong_vtu, vertical, 1)
_, *weak_v = sample(args.weak_vtu, vertical, 1)
_, *nosupg_v = sample(args.weak_nosupg_vtu, vertical, 1)
write_dat(out_dir / f"{mesh_stem}_{re_label}_vertical.dat", "y", y, strong_v, weak_v, nosupg_v)
