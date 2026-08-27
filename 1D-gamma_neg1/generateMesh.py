import gmsh
import argparse

parser = argparse.ArgumentParser(
    prog="mesh generator",
    description="Creates a 1D mesh with nx divisions",
)

parser.add_argument("-nx", "--numberElements", type=int)
args = parser.parse_args()
nx = args.numberElements
lc = 0.1

gmsh.initialize()

gmsh.model.add("domain")

gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(1, 0, 0, lc, 2)
gmsh.model.geo.addLine(1, 2, 1)

gmsh.model.geo.mesh.setTransfiniteCurve(1, nx + 1)
gmsh.model.geo.synchronize()

gmsh.model.addPhysicalGroup(0, [1], name="inlet")
gmsh.model.addPhysicalGroup(0, [2], name="outlet")
gmsh.model.addPhysicalGroup(1, [1], name="domain")

gmsh.model.mesh.generate(1)

gmsh.write(f"meshes/{nx}-line.msh")

gmsh.finalize()
