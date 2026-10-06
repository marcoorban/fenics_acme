from pathlib import Path

from dolfinx.io import VTXWriter


class VTX_Files:
    def __init__(self, folder):
        self.folder = Path(folder)
        self.vtx_u = None
        self.vtx_p = None

    def __call__(self, t, u_, p_):
        if self.vtx_u is None:
            comm = u_.function_space.mesh.comm
            self.vtx_u = VTXWriter(comm, self.folder / "u.bp", [u_], engine="BP4")
            self.vtx_p = VTXWriter(comm, self.folder / "p.bp", [p_], engine="BP4")
        self.vtx_u.write(t)
        self.vtx_p.write(t)

    def close(self):
        if self.vtx_u is not None:
            self.vtx_u.close()
            self.vtx_p.close()
