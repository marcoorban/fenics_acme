#!/usr/bin/env python3
"""
Generate a structured (tensor-product) quadrilateral mesh of a rectangle,
written directly as a Gmsh 4.1 ASCII .msh file.

The geometry and physical groups reproduce those of the original mesh:

      (0,1) 4 ------- 3 (1,1)
            |    3    |
       (4)  |         | (2)      curve 1 -> physical 2  "bottom"
   y_split 5|  domain |          curve 2 -> physical 3  "right"
       (5)  |   (10)  |          curve 3 -> physical 4  "top"
            |         |          curve 4 -> physical 5  "left_top"
      (0,0) 1 ------- 2 (1,0)    curve 5 -> physical 1  "left_bottom"
                 1                surface 1 -> physical 10 "domain"

The left edge is split at y = y_split, so that node line must exist in the
grid: ny * y_split must be an integer.

Usage:
    python structured_mesh.py [-o out.msh] [--nx NX] [--ny NY]
                              [--lx LX] [--ly LY] [--y-split YS]
"""

import argparse


def build_mesh(nx, ny, lx=1.0, ly=1.0, y_split=0.4):
    """Return (nodes, quads, edges) for a structured nx-by-ny grid.

    nodes : dict (i, j) -> global 1-based tag
    """
    js_f = ny * (y_split / ly)
    js = round(js_f)
    if abs(js_f - js) > 1e-12:
        raise ValueError(
            f"The split point y={y_split} does not fall on a grid line for "
            f"ny={ny} (ny*y_split/ly = {js_f}). Choose ny as a multiple of "
            f"{int(round(1 / (y_split / ly)))} or adjust y_split."
        )
    if not 0 < js < ny:
        raise ValueError("y_split must lie strictly inside (0, ly).")

    tag = {}  # (i, j) -> node tag
    blocks = []  # (dim, entity_tag, [(i, j), ...]) in node-output order
    counter = 0

    def add(entity_dim, entity_tag, ij_list):
        nonlocal counter
        for ij in ij_list:
            counter += 1
            tag[ij] = counter
        blocks.append((entity_dim, entity_tag, ij_list))

    # --- dim 0: the five geometric points (corners + split point) ---------
    add(0, 1, [(0, 0)])
    add(0, 2, [(nx, 0)])
    add(0, 3, [(nx, ny)])
    add(0, 4, [(0, ny)])
    add(0, 5, [(0, js)])

    # --- dim 1: interior nodes of each curve, in curve direction ----------
    add(1, 1, [(i, 0) for i in range(1, nx)])  # bottom  1->2
    add(1, 2, [(nx, j) for j in range(1, ny)])  # right   2->3
    add(1, 3, [(i, ny) for i in range(nx - 1, 0, -1)])  # top     3->4
    add(1, 4, [(0, j) for j in range(ny - 1, js, -1)])  # l_top   4->5
    add(1, 5, [(0, j) for j in range(js - 1, 0, -1)])  # l_bot   5->1

    # --- dim 2: interior nodes of the surface -----------------------------
    add(2, 1, [(i, j) for j in range(1, ny) for i in range(1, nx)])

    assert counter == (nx + 1) * (ny + 1)

    # --- cells: quads, counter-clockwise ----------------------------------
    quads = [
        (tag[(i, j)], tag[(i + 1, j)], tag[(i + 1, j + 1)], tag[(i, j + 1)])
        for j in range(ny)
        for i in range(nx)
    ]

    # --- boundary line elements, one block per curve ----------------------
    def path(ij_seq):
        return [(tag[a], tag[b]) for a, b in zip(ij_seq, ij_seq[1:])]

    edges = [
        (1, path([(i, 0) for i in range(nx + 1)])),
        (2, path([(nx, j) for j in range(ny + 1)])),
        (3, path([(i, ny) for i in range(nx, -1, -1)])),
        (4, path([(0, j) for j in range(ny, js - 1, -1)])),
        (5, path([(0, j) for j in range(js, -1, -1)])),
    ]

    return blocks, tag, quads, edges, js


def write_msh(path, nx, ny, lx=1.0, ly=1.0, y_split=0.4):
    blocks, tag, quads, edges, js = build_mesh(nx, ny, lx, ly, y_split)
    n_nodes = (nx + 1) * (ny + 1)
    n_elems = len(quads) + sum(len(e) for _, e in edges)

    def fmt(v):
        return repr(float(v))

    def xy(ij):
        i, j = ij
        return fmt(lx * i / nx), fmt(ly * j / ny)

    out = []
    w = out.append

    w("$MeshFormat\n4.1 0 8\n$EndMeshFormat")

    w("$PhysicalNames\n6")
    for dim, ptag, name in [
        (1, 1, "left_bottom"),
        (1, 2, "bottom"),
        (1, 3, "right"),
        (1, 4, "top"),
        (1, 5, "left_top"),
        (2, 10, "domain"),
    ]:
        w(f'{dim} {ptag} "{name}"')
    w("$EndPhysicalNames")

    # --- $Entities: 5 points, 5 curves, 1 surface, 0 volumes --------------
    ys = fmt(ly * js / ny)
    LX, LY, Z = fmt(lx), fmt(ly), "0"
    w("$Entities\n5 5 1 0")
    w(f"1 0 0 0 0")  # point tags carry no physical group
    w(f"2 {LX} 0 0 0")
    w(f"3 {LX} {LY} 0 0")
    w(f"4 0 {LY} 0 0")
    w(f"5 0 {ys} 0 0")
    #    tag  minX minY minZ  maxX maxY maxZ  nPhys phys  nBndPts  pts
    w(f"1 0 0 0 {LX} 0 0 1 2 2 1 -2")
    w(f"2 {LX} 0 0 {LX} {LY} 0 1 3 2 2 -3")
    w(f"3 0 {LY} 0 {LX} {LY} 0 1 4 2 3 -4")
    w(f"4 0 {ys} 0 0 {LY} 0 1 5 2 4 -5")
    w(f"5 0 0 0 0 {ys} 0 1 1 2 5 -1")
    w(f"1 0 0 0 {LX} {LY} 0 1 10 5 1 2 3 4 5")
    w("$EndEntities")

    # --- $Nodes -----------------------------------------------------------
    w(f"$Nodes\n{len(blocks)} {n_nodes} 1 {n_nodes}")
    for dim, etag, ijs in blocks:
        w(f"{dim} {etag} 0 {len(ijs)}")
        for ij in ijs:
            w(str(tag[ij]))
        for ij in ijs:
            x, y = xy(ij)
            w(f"{x} {y} 0")
    w("$EndNodes")

    # --- $Elements --------------------------------------------------------
    w(f"$Elements\n{len(edges) + 1} {n_elems} 1 {n_elems}")
    eid = 0
    for ctag, segs in edges:  # type 1 = 2-node line
        w(f"1 {ctag} 1 {len(segs)}")
        for a, b in segs:
            eid += 1
            w(f"{eid} {a} {b}")
    w(f"2 1 3 {len(quads)}")  # type 3 = 4-node quad
    for q in quads:
        eid += 1
        w(f"{eid} {q[0]} {q[1]} {q[2]} {q[3]}")
    w("$EndElements")

    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")

    return n_nodes, len(quads), n_elems


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("-o", "--output", default="rectangle.msh")
    p.add_argument("--nx", type=int, default=65, help="cells in x (default 65)")
    p.add_argument("--ny", type=int, default=65, help="cells in y (default 65)")
    p.add_argument("--lx", type=float, default=1.0)
    p.add_argument("--ly", type=float, default=1.0)
    p.add_argument(
        "--y-split",
        type=float,
        default=0.4,
        help="left boundary split height (default 0.4)",
    )
    a = p.parse_args()

    outputName = f"rectangle-{a.nx}x{a.ny}.msh"
    nn, nq, ne = write_msh(outputName, a.nx, a.ny, a.lx, a.ly, a.y_split)
    print(f"{a.output}: {a.nx}x{a.ny} structured quad grid")
    print(f"  nodes {nn}, quads {nq}, elements {ne} (incl. boundary lines)")
    print(f"  hx = {a.lx / a.nx:.6g}, hy = {a.ly / a.ny:.6g}")


if __name__ == "__main__":
    main()
