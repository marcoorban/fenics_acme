import argparse

import gmsh

parser = argparse.ArgumentParser(
    description="Generate a structured NxN-element cavity mesh for the FEM "
    "solver -- same lid/wall/domain physical-group convention as "
    "navStokes/cavity/mesh/mesh.py (tags 1/2/3), without that script's "
    "OpenFOAM 3D extrusion. Deliberately coarse diagnostic mesh, not a "
    "replacement for mesh.py."
)
parser.add_argument("--n", type=int, default=10, help="Elements per side.")
parser.add_argument(
    "--bump",
    type=float,
    default=1.0,
    help="Transfinite 'Bump' coefficient (ratio of largest to smallest cell "
    "on each curve). 1.0 = uniform spacing (default). Smaller = stronger "
    "clustering toward the two ends of each curve, i.e. toward all four "
    "walls (same convention as mesh.py's bumpCoef=0.15), at the same total "
    "element count -- redistributing resolution rather than adding more.",
)
parser.add_argument("--out", default=None, help="Output .msh path.")
args = parser.parse_args()
out = args.out or (
    f"cavity_structured_{args.n}x{args.n}.msh"
    if args.bump == 1.0
    else f"cavity_structured_{args.n}x{args.n}_bump{args.bump:g}.msh"
)

d = 1
gmsh.initialize()

lowerLeft = gmsh.model.geo.addPoint(0, 0, 0)
lowerRight = gmsh.model.geo.addPoint(d, 0, 0)
upperRight = gmsh.model.geo.addPoint(d, d, 0)
upperLeft = gmsh.model.geo.addPoint(0, d, 0)
bottom = gmsh.model.geo.addLine(lowerLeft, lowerRight)
right = gmsh.model.geo.addLine(lowerRight, upperRight)
top = gmsh.model.geo.addLine(upperRight, upperLeft)
left = gmsh.model.geo.addLine(upperLeft, lowerLeft)
loop = gmsh.model.geo.addCurveLoop([bottom, right, top, left])
surface = gmsh.model.geo.addPlaneSurface([loop])

# n elements per side means n+1 points per side. "Bump" with coef=1.0
# degenerates to uniform spacing, same as mesh.py's own comment notes.
points = args.n + 1
for curve in (bottom, top, left, right):
    gmsh.model.geo.mesh.setTransfiniteCurve(
        curve, points, meshType="Bump", coef=args.bump
    )
gmsh.model.geo.mesh.setTransfiniteSurface(surface)
gmsh.model.geo.mesh.setRecombine(2, surface)

gmsh.model.geo.synchronize()

# Same creation order as mesh.py -> same auto-assigned tags (1=lid, 2=wall,
# 3=domain), which fem_*.py's ft.find(1)/ft.find(2) rely on.
gmsh.model.addPhysicalGroup(1, [top], name="lid")
gmsh.model.addPhysicalGroup(1, [bottom, right, left], name="wall")
gmsh.model.addPhysicalGroup(2, [surface], name="domain")

gmsh.model.mesh.generate(2)
gmsh.write(out)
gmsh.finalize()

print(f"Wrote {out} ({args.n}x{args.n} = {args.n * args.n} elements)")
