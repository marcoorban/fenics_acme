import argparse
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

parser = argparse.ArgumentParser(
    description="Check mass conservation (continuity, div(u)=0) per mesh "
    "element from a saved velocity .vtu file."
)
parser.add_argument("u_vtu", help="Path to the velocity .vtu file (e.g. results/u_p0_000000.vtu).")
args = parser.parse_args()

reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName(args.u_vtu)
reader.Update()
grid = reader.GetOutput()

# Gradient of the velocity field, then its trace (divergence), computed
# from the actual (quadratic-Lagrange) point data -- not re-solving
# anything, just differentiating what was written to disk.
gradient = vtk.vtkGradientFilter()
gradient.SetInputData(grid)
gradient.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "u")
gradient.SetComputeDivergence(True)
gradient.SetDivergenceArrayName("div_u")
gradient.Update()

# Average the point-wise divergence onto one value per element.
point_to_cell = vtk.vtkPointDataToCellData()
point_to_cell.SetInputData(gradient.GetOutput())
point_to_cell.Update()
out = point_to_cell.GetOutput()

div_u = vtk_to_numpy(out.GetCellData().GetArray("div_u"))
print(f"Cells: {len(div_u)}")
print(f"max|div(u)|  = {np.max(np.abs(div_u))}")
print(f"mean|div(u)| = {np.mean(np.abs(div_u))}")
print(f"RMS div(u)   = {np.sqrt(np.mean(div_u**2))}")

out_path = Path(args.u_vtu).with_name("continuity_check.vtu")
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName(str(out_path))
writer.SetInputData(out)
writer.Write()
print(f"Wrote per-element divergence field -> {out_path}")
