from pathlib import Path

from dolfinx.io import gmsh
from mpi4py import MPI


def read_mesh(fileName):
    cylinderDir = Path(__file__).parent.parent.parent
    meshFile = cylinderDir / "mesh" / fileName
    meshData = gmsh.read_from_msh(meshFile, MPI.COMM_WORLD, rank=0, gdim=2)
    mesh = meshData.mesh
    assert meshData.facet_tags is not None
    ft = meshData.facet_tags
    ft.name = "Facet markers"
    return (mesh, ft)



