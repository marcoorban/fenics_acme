import gmsh
from mpi4py import MPI

gmsh.initialize()

comm = MPI.COMM_WORLD
model_rank = 0

if comm.rank == model_rank:
    gmsh.model.add("rectangle")
    geo = gmsh.model.geo

    L, H = 1.0, 1.0
    Y = 0.4
    nx = 65
    ny = 65
    ny_bottom = int(Y * 65)
    ny_top = ny - ny_bottom

    p1 = geo.addPoint(0.0, 0.0, 0.0)
    p2 = geo.addPoint(L, 0.0, 0.0)
    p3 = geo.addPoint(L, H, 0.0)
    p4 = geo.addPoint(0.0, H, 0.0)
    p5 = geo.addPoint(0.0, Y, 0.0)

    l_bottom = geo.addLine(p1, p2)
    l_right = geo.addLine(p2, p3)
    l_top = geo.addLine(p3, p4)
    l_left_top = geo.addLine(p4, p5)
    l_left_bottom = geo.addLine(p5, p1)

    loop = geo.addCurveLoop([l_bottom, l_right, l_top, l_left_top, l_left_bottom])
    surf = geo.addPlaneSurface([loop])

    # transfite constraints  #
    geo.mesh.setTransfiniteCurve(l_bottom, nx)
    geo.mesh.setTransfiniteCurve(l_top, nx)
    geo.mesh.setTransfiniteCurve(l_left_top, ny_top)
    geo.mesh.setTransfiniteCurve(l_left_bottom, ny_bottom)
    geo.mesh.setTransfiniteCurve(l_right, ny)

    geo.mesh.setTransfiniteSurface(surf, "Left", [p1, p2, p3, p4])

    geo.mesh.setRecombine(2, surf)

    geo.synchronize()

gmsh.model.addPhysicalGroup(1, [l_left_bottom], tag=1, name="left_bottom")
gmsh.model.addPhysicalGroup(1, [l_bottom], tag=2, name="bottom")
gmsh.model.addPhysicalGroup(1, [l_right], tag=3, name="right")
gmsh.model.addPhysicalGroup(1, [l_top], tag=4, name="top")
gmsh.model.addPhysicalGroup(1, [l_left_top], tag=5, name="left_top")

# Tag the surface also
gmsh.model.addPhysicalGroup(2, [surf], tag=10, name="domain")

## Step 3 - Mesh
gmsh.model.mesh.generate(2)
gmsh.model.mesh.setOrder(1)
gmsh.write("rectangle.msh")
gmsh.finalize()
