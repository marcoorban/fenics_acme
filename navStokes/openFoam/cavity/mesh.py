import gmsh 
import sys

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

# Add transfinite points
points = 50
gmsh.model.geo.mesh.setTransfiniteCurve(bottom, points)
gmsh.model.geo.mesh.setTransfiniteCurve(top, points)
gmsh.model.geo.mesh.setTransfiniteCurve(left, points)
gmsh.model.geo.mesh.setTransfiniteCurve(right, points)
gmsh.model.geo.mesh.setTransfiniteSurface(surface)

gmsh.model.geo.mesh.setRecombine(2, surface)

# Extrude one unit in z to give OpenFoam a 3D mesh.
# Extrude returns a list with the created entities in the following order:
#   out[0] = far cap (z=d), out[1] = volume,
#   out[2..5] = side surfaces extruded from bottom, right, top, left
out = gmsh.model.geo.extrude([(2, surface)], 0, 0, d, numElements=[1], recombine=True)
farCap, volume, bottomWall, rightWall, lidWall, leftWall = out

gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(2, [lidWall[1]], name="lid")
gmsh.model.addPhysicalGroup(2, [bottomWall[1], rightWall[1], leftWall[1]], name="wall")
# Front/back caps carry no physics in a 2D flow problem.
gmsh.model.addPhysicalGroup(2, [surface, farCap[1]], name="empty")
gmsh.model.addPhysicalGroup(3, [volume[1]], name="domain")

gmsh.model.mesh.generate(3)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("cavity.msh")
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
gmsh.finalize()

