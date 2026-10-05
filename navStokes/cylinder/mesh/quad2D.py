import gmsh
import sys

gmsh.initialize()

# --- Parameters ---
D = 0.1  # cylinder diameter
R = D / 2


thickness = 0.1 * D  # thin single-layer extrusion for OpenFOAM's pseudo-2D
# convention (see cavity/mesh/mesh.py) -- not physically meaningful, just
# gives OpenFOAM the 3D mesh it requires even for a 2D problem.

# FEM (Re=100, dokken/navStokes.py) mesh sizing -- Distance+Threshold field
# grades the size smoothly from size_cyl at the cylinder surface out to
# size_far beyond dist_max. Coarsened 2026-09-30 to the dokken tutorial's own
# resolution (~D/6 to D/3 at the cylinder, i.e. ~0.0167 for D=0.1) after the
# previous values -- leftover from the Re=20 OpenFOAM 2D-1 sizing study,
# h_min~0.0016 once meshed -- caused a Courant-number blow-up with the
# tutorial's dt=1/1600 (see cylinder/log.md, 2026-09-30 #fem entry).
size_far = 0.08  # target element size far from the cylinder -- bumped
# 2026-09-30 (was 0.017, still visibly too fine near the domain ends/outlet);
# closer to the dokken tutorial's own far-field size (LcMax=0.25*H~0.1025).
size_cyl = 0.015  # target element size right at the cylinder surface
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

# Rectangle with the cylinder cut out as a hole -- unstructured quads
# (recombined below), refined near the cylinder via the Distance+Threshold
# field below.
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

# Quadrilateral elements -- 2026-09-30: the dokken tutorial's own mesh (see
# navstokes_original.py) uses quads and runs perfectly at dt=1/1600 on this
# same problem; this triangular mesh (even after coarsening to the same
# element size) still blew up, isolating the element type/family -- not just
# the size -- as the actual problem. Same settings the tutorial uses:
# Frontal-Delaunay-for-Quads meshing algorithm + full recombination into
# quads, applied on top of the existing graded Distance+Threshold sizing
# field above (fine at the cylinder, coarse far away is preserved either way
# -- recombination pairs up the already-graded triangles, it doesn't re-size
# them).
gmsh.option.setNumber("Mesh.Algorithm", 8)
gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
gmsh.model.geo.mesh.setRecombine(2, surface)

# Extrude one thin layer in z. Extrude returns entities in the order:
#   out[0] = far cap, out[1] = volume,
#   out[2..] = lateral surfaces in curve-loop order -- outer loop first
#   (bottom, outlet, top, inlet), then the cylinder loop (its 4 arcs).
# out = gmsh.model.geo.extrude(
#    [(2, surface)], 0, 0, thickness, numElements=[1], recombine=True
# )
# farCap, volume = out[0], out[1]
# bottomWall, outletFace, topWall, inletFace = out[2], out[3], out[4], out[5]
# cylArcs = [out[6][1], out[7][1], out[8][1], out[9][1]]

# gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(1, [inlet], 1, name="inlet")
gmsh.model.addPhysicalGroup(1, [outlet], 2, name="outlet")
gmsh.model.addPhysicalGroup(1, [top, bottom], 3, name="wall")
# 2026-09-30: was `[cylinderLoop]` -- a curve-LOOP tag (from addCurveLoop),
# not a curve entity; loops aren't valid members of a dim=1 physical group,
# so that silently produced an EMPTY group (`ft.find(4)` -> 0 facets). The
# no-slip cylinder BC in navStokes.py was therefore never actually applied --
# the hole in the domain was a real geometric hole, but its boundary was an
# unconstrained "do-nothing" boundary instead of a wall. Present since this
# file's first commit; likely the real cause of the persistent blow-up
# (mesh-size/Courant fixes alone couldn't have fixed a missing wall BC).
gmsh.model.addPhysicalGroup(
    1, [arcRightTop, arcTopLeft, arcLeftBottom, arcBottomRight], 4, name="cylinder"
)
gmsh.model.addPhysicalGroup(2, [surface], 5, name="domain")

# generate(2), not generate(3) -- the extrude block above is commented out,
# so there are no 3D entities; this mesh is genuinely 2D, matching the
# "domain" physical group above (dim=2, the surface itself).
gmsh.model.mesh.generate(2)

# Filename is parametrized by the actual generated cell count, not hardcoded
# -- read from gmsh's own result (not size_cyl/size_far, which only predict
# it approximately, see cylinder/log.md's mesh-sizing-iteration entry) so it
# can't drift stale if the sizing parameters above are retuned. This also
# doubles as this mesh's id in cylinder/mesh/meshes_summary.csv -- re-run
# that catalog's numbers if this changes.
elem_types, elem_tags, _ = gmsh.model.mesh.getElements(dim=2)
num_cells = sum(len(tags) for tags in elem_tags)
meshFileName = f"quad2D_{num_cells}cells.msh"
print(f"Generated {num_cells} cells -> writing {meshFileName}")
gmsh.write(meshFileName)
if "-nopopup" not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()
