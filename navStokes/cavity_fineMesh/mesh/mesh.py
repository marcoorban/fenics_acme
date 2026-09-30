import gmsh 
import sys
from pathlib import Path 
import os

scriptDir = os.path.dirname(os.path.realpath(__file__))
fileName = "cavity.msh"
meshFile = str(scriptDir) + "/" + fileName
gmsh.initialize()

d = 1
lowerLeft = gmsh.model.geo.addPoint(0, 0, 0)
lowerRight = gmsh.model.geo.addPoint(d, 0, 0)
upperRight = gmsh.model.geo.addPoint(d, d, 0)
upperLeft = gmsh.model.geo.addPoint(0, d, 0)
# Connect them with lines 
bottom = gmsh.model.geo.addLine(lowerLeft, lowerRight, 11)
right = gmsh.model.geo.addLine(lowerRight, upperRight, 12)
top = gmsh.model.geo.addLine(upperRight, upperLeft, 13)
left = gmsh.model.geo.addLine(upperLeft, lowerLeft, 14)
# Create a loop
loop = gmsh.model.geo.addCurveLoop([bottom, right, top, left], 15)
# Create a surface
surface = gmsh.model.geo.addPlaneSurface([loop], 21)

# Corner-vortex refinement, replacing the old uniform transfinite mesh with
# a Distance+Threshold field (same pattern as cylinder/mesh/mesh.py). The
# classic lid-driven-cavity corner (Moffatt) eddies are geometrically
# self-similar and shrink rapidly toward each corner -- resolving them needs
# cells far smaller than a *uniform* mesh fine enough for the bulk flow
# would ever need, and a uniform mesh at corner resolution everywhere would
# be wildly expensive. Grade smoothly instead: size_wall right at all four
# walls/corners (where corner eddies live -- not just the bottom corners,
# since secondary eddies can appear near the top corners too at higher Re),
# size_far in the bulk interior, linear transition between dist_min/dist_max.
# First-pass values, not yet validated against an actual solve -- iterate
# empirically (checkMesh cell count / visual eddy resolution) like the
# cylinder mesh's sizing was tuned.
size_wall = 0.002  # target element size at the cavity walls/corners
size_far = 0.02  # target element size in the bulk interior
dist_min = 0.05  # distance from a wall where size_wall still fully applies
dist_max = 0.3  # distance beyond which size_far fully applies (linear
# grading from size_wall to size_far between dist_min and dist_max)

# Curves must be synchronized into the model before a field can reference
# them; extrude() below adds more entities and is synchronized separately.
gmsh.model.geo.synchronize()

distField = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(distField, "CurvesList", [bottom, right, top, left])
gmsh.model.mesh.field.setNumber(distField, "Sampling", 200)

thresholdField = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(thresholdField, "InField", distField)
gmsh.model.mesh.field.setNumber(thresholdField, "SizeMin", size_wall)
gmsh.model.mesh.field.setNumber(thresholdField, "SizeMax", size_far)
gmsh.model.mesh.field.setNumber(thresholdField, "DistMin", dist_min)
gmsh.model.mesh.field.setNumber(thresholdField, "DistMax", dist_max)

gmsh.model.mesh.field.setAsBackgroundMesh(thresholdField)

# Let the field be the only sizing source -- otherwise gmsh blends it with
# point-/curvature-based sizing, which can reintroduce sharp transitions
# (same lesson as cylinder/mesh/mesh.py).
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

# The free (Frontal-Delaunay) algorithm triangulates by default; recombine
# the base surface into quads to keep the structured-like cell shape the
# solver was tuned against, same as the transfinite mesh gave before.
gmsh.model.geo.mesh.setRecombine(2, surface)

# Extrude one unit in z to give OpenFoam a 3D mesh.
# Extrude returns a list with the created entities in the following order:
#   out[0] = far cap (z=d), out[1] = volume,
#   out[2..5] = side surfaces extruded from bottom, right, top, left
out = gmsh.model.geo.extrude([(2, surface)], 0, 0, d/100, numElements=[1], recombine=True)
farCap, volume, bottomWall, rightWall, lidWall, leftWall = out

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(2, [lidWall[1]], name="lid")
gmsh.model.addPhysicalGroup(2, [bottomWall[1], rightWall[1], leftWall[1]], name="wall")
# Front/back caps carry no physics in a 2D flow problem.
gmsh.model.addPhysicalGroup(2, [surface, farCap[1]], name="empty")
gmsh.model.addPhysicalGroup(3, [volume[1]], name="domain")

gmsh.model.mesh.generate(3)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write(meshFile)
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()

