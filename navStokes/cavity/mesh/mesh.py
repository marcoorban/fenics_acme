import gmsh 
import sys
from pathlib import Path 
import os

scriptDir = os.path.dirname(os.path.realpath(__file__))
fileName = "cavity.msh"
meshFile = str(scriptDir) + "/" + fileName
fileName2D = "cavity-2D.msh"
meshFile2D = str(scriptDir) + "/" + fileName2D
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

# Add transfinite points, clustered toward both ends of each curve ("Bump"
# distribution) instead of uniform spacing -- since every wall is an
# endpoint of two of these four curves (e.g. the bottom-left corner is the
# shared endpoint of "bottom" and "left"), biasing all four curves this way
# refines the grid near all four cavity walls/corners and leaves the bulk
# interior coarser, without changing the point count (still 130 points,
# i.e. 129 cells, per edge -- same total cell count as the uniform grid).
# coef is the ratio of the largest to smallest cell on the curve; smaller
# coef => stronger clustering at the two ends.
points = 130
bumpCoef = 0.15
gmsh.model.geo.mesh.setTransfiniteCurve(bottom, points, meshType="Bump", coef=bumpCoef)
gmsh.model.geo.mesh.setTransfiniteCurve(top, points, meshType="Bump", coef=bumpCoef)
gmsh.model.geo.mesh.setTransfiniteCurve(left, points, meshType="Bump", coef=bumpCoef)
gmsh.model.geo.mesh.setTransfiniteCurve(right, points, meshType="Bump", coef=bumpCoef)
gmsh.model.geo.mesh.setTransfiniteSurface(surface)

gmsh.model.geo.mesh.setRecombine(2, surface)

# Save the plain 2D mesh (no extrusion) for FEniCS, which needs a genuinely
# 2D mesh -- unlike OpenFOAM, it has no use for the pseudo-2D single-layer
# extrusion below. Must synchronize + generate(2) here, before extruding:
# once the surface is extruded into a volume, generate(3) meshes it as part
# of a 3D geometry and there's no clean "2D-only" mesh left to write out.
gmsh.model.geo.synchronize()

# Same physical-group convention as the 3D mesh below, minus "empty" -- that
# group only marks OpenFOAM's pseudo-2D front/back caps as non-physical,
# which don't exist on a genuinely 2D mesh. "domain" here is the surface
# itself (the top-dimensional entity), not a volume.
gmsh.model.addPhysicalGroup(1, [top], name="lid")
gmsh.model.addPhysicalGroup(1, [bottom, right, left], name="wall")
gmsh.model.addPhysicalGroup(2, [surface], name="domain")

gmsh.model.mesh.generate(2)
# Written with gmsh's current default (latest) MSH format -- only the
# OpenFOAM mesh below needs the old version 2.2 format, set just before
# that write, so it must not be set yet at this point.
gmsh.write(meshFile2D)

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
gmsh.finalize()

