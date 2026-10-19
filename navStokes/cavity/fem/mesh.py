import gmsh


def create_mesh():
    """Unit square [0,1]x[0,1], quad mesh, refined near all four edges
    (graded via a Distance-to-boundary + Threshold sizing field). Writes
    the mesh to a .msh file instead of returning an in-memory dolfinx
    mesh -- run this once, then have the solver read the saved file
    (see meshreader.py) rather than regenerating the mesh every run."""
    gmsh.initialize()
    gmsh.model.add("cavity")

    p1 = gmsh.model.geo.addPoint(0, 0, 0)
    p2 = gmsh.model.geo.addPoint(1, 0, 0)
    p3 = gmsh.model.geo.addPoint(1, 1, 0)
    p4 = gmsh.model.geo.addPoint(0, 1, 0)

    bottom = gmsh.model.geo.addLine(p1, p2)
    right = gmsh.model.geo.addLine(p2, p3)
    top = gmsh.model.geo.addLine(p3, p4)
    left = gmsh.model.geo.addLine(p4, p1)

    loop = gmsh.model.geo.addCurveLoop([bottom, right, top, left])
    surface = gmsh.model.geo.addPlaneSurface([loop])
    gmsh.model.geo.synchronize()

    # Grade element size up from the four edges toward the center.
    dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(dist, "CurvesList", [bottom, right, top, left])

    size = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(size, "InField", dist)
    gmsh.model.mesh.field.setNumber(size, "SizeMin", 0.01)
    gmsh.model.mesh.field.setNumber(size, "SizeMax", 0.05)
    gmsh.model.mesh.field.setNumber(size, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(size, "DistMax", 0.3)
    gmsh.model.mesh.field.setAsBackgroundMesh(size)

    # Quadrilateral elements: mesh with triangles first (Frontal-Delaunay-
    # for-Quads), then recombine pairs of triangles into quads. The sizing
    # field above still controls density -- recombination pairs up
    # already-graded triangles, it doesn't re-size them.
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 2)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 1)
    gmsh.model.geo.mesh.setRecombine(2, surface)

    # Top edge is the moving lid; the other three are the fixed walls.
    gmsh.model.addPhysicalGroup(1, [top], 1, name="lid")
    gmsh.model.addPhysicalGroup(1, [bottom, left, right], 2, name="wall")
    gmsh.model.addPhysicalGroup(2, [surface], 3, name="domain")

    gmsh.model.mesh.generate(2)

    # Filename carries the actual generated cell count, read back from
    # gmsh's own result (same convention as cylinder/mesh/quad2D.py).
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements(dim=2)
    num_cells = sum(len(tags) for tags in elem_tags)
    meshFileName = f"cavity_{num_cells}cells.msh"
    gmsh.write(meshFileName)
    print(f"Generated {num_cells} cells -> wrote {meshFileName}")
    gmsh.finalize()


if __name__ == "__main__":
    create_mesh()
