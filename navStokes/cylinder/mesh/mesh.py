import gmsh
import sys

gmsh.initialize()

# --- Parameters ---
D = 0.1  # cylinder diameter
R = D / 2


thickness = 0.1 * D  # thin single-layer extrusion for OpenFOAM's pseudo-2D
# convention (see cavity/mesh/mesh.py) -- not physically meaningful, just
# gives OpenFOAM the 3D mesh it requires even for a 2D problem.

# Deliberately coarse for now -- the previous per-point sizing (mesh_size_far
# at every far point, mesh_size_cyl right at the cylinder, nothing in
# between) jumped straight from 0.01 to 0.001 with no graded transition,
# which checkMesh flagged as ~12% severely non-orthogonal (>70 deg) faces
# and which correlated with a Courant-number blow-up during solving (see
# cylinder/log.md, 2026-09-10). Below, a Distance+Threshold field grades
# the size smoothly instead of jumping. These values aren't meant to give
# an accurate solution -- just a mesh coarse/well-behaved enough to run to
# completion without blowing up; refine later once that's confirmed.
size_far = 0.0097  # target element size far from the cylinder
size_cyl = 0.00194  # target element size right at the cylinder surface
dist_min = D  # distance from the cylinder surface where size_cyl still applies
dist_max = 3 * D  # distance beyond which size_far fully applies (linear
# grading from size_cyl to size_far between dist_min and dist_max)

# --- Outer rectangle ---
x_min = 0
x_max = 2.2
y_min = 0
y_max = 0.41

inlet_dist = 0.15
vertical_offset = 0.15
cx, cy = inlet_dist + R, vertical_offset + R

p1 = gmsh.model.geo.addPoint(x_min, y_min, 0)
p2 = gmsh.model.geo.addPoint(x_max, y_min, 0)
p3 = gmsh.model.geo.addPoint(x_max, y_max, 0)
p4 = gmsh.model.geo.addPoint(x_min, y_max, 0)

bottom = gmsh.model.geo.addLine(p1, p2)
outlet = gmsh.model.geo.addLine(p2, p3)
top = gmsh.model.geo.addLine(p3, p4)
inlet = gmsh.model.geo.addLine(p4, p1)

outerLoop = gmsh.model.geo.addCurveLoop([bottom, outlet, top, inlet])

# --- Cylinder ---
pc = gmsh.model.geo.addPoint(cx, cy, 0)
pRight = gmsh.model.geo.addPoint(cx + R, cy, 0)
pTop = gmsh.model.geo.addPoint(cx, cy + R, 0)
pLeft = gmsh.model.geo.addPoint(cx - R, cy, 0)
pBottom = gmsh.model.geo.addPoint(cx, cy - R, 0)

arcRightTop = gmsh.model.geo.addCircleArc(pRight, pc, pTop)
arcTopLeft = gmsh.model.geo.addCircleArc(pTop, pc, pLeft)
arcLeftBottom = gmsh.model.geo.addCircleArc(pLeft, pc, pBottom)
arcBottomRight = gmsh.model.geo.addCircleArc(pBottom, pc, pRight)

cylinderLoop = gmsh.model.geo.addCurveLoop(
    [arcRightTop, arcTopLeft, arcLeftBottom, arcBottomRight]
)

# Rectangle with the cylinder cut out as a hole -- unstructured triangles,
# refined near the cylinder via the per-point mesh sizes above.
surface = gmsh.model.geo.addPlaneSurface([outerLoop, cylinderLoop])

gmsh.model.geo.synchronize()

# Grade element size smoothly from size_cyl (at the cylinder surface) out to
# size_far (beyond dist_max), instead of the old sharp per-point jump.
distField = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(
    distField,
    "CurvesList",
    [arcRightTop, arcTopLeft, arcLeftBottom, arcBottomRight],
)
gmsh.model.mesh.field.setNumber(distField, "Sampling", 100)

thresholdField = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(thresholdField, "InField", distField)
gmsh.model.mesh.field.setNumber(thresholdField, "SizeMin", size_cyl)
gmsh.model.mesh.field.setNumber(thresholdField, "SizeMax", size_far)
gmsh.model.mesh.field.setNumber(thresholdField, "DistMin", dist_min)
gmsh.model.mesh.field.setNumber(thresholdField, "DistMax", dist_max)

gmsh.model.mesh.field.setAsBackgroundMesh(thresholdField)

# Let the field be the only source of sizing -- otherwise gmsh blends it
# with curvature-based/point-based sizing in ways that can reintroduce the
# sharp transitions this field is meant to avoid.
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

# Extrude one thin layer in z. Extrude returns entities in the order:
#   out[0] = far cap, out[1] = volume,
#   out[2..] = lateral surfaces in curve-loop order -- outer loop first
#   (bottom, outlet, top, inlet), then the cylinder loop (its 4 arcs).
out = gmsh.model.geo.extrude(
    [(2, surface)], 0, 0, thickness, numElements=[1], recombine=True
)
farCap, volume = out[0], out[1]
bottomWall, outletFace, topWall, inletFace = out[2], out[3], out[4], out[5]
cylArcs = [out[6][1], out[7][1], out[8][1], out[9][1]]

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(2, [inletFace[1]], name="inlet")
gmsh.model.addPhysicalGroup(2, [outletFace[1]], name="outlet")
gmsh.model.addPhysicalGroup(2, [topWall[1], bottomWall[1]], name="wall")
gmsh.model.addPhysicalGroup(2, cylArcs, name="cylinder")
gmsh.model.addPhysicalGroup(2, [surface, farCap[1]], name="empty")
gmsh.model.addPhysicalGroup(3, [volume[1]], name="domain")

gmsh.model.mesh.generate(3)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("cylinder.msh")
if "-nopopup" not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
