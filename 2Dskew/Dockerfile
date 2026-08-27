FROM dolfinx/dolfinx:nightly

# gmsh's Python wheel links against GUI libs it never uses headlessly;
# without them `import gmsh` fails with "libGLU.so.1: cannot open shared
# object file" even though dolfinx.io.gmsh only reads a mesh file with it.
RUN apt-get update && apt-get install -y --no-install-recommends \
    libglu1-mesa \
    libxrender1 \
    libxcursor1 \
    libxinerama1 \
    libxft2 \
    && rm -rf /var/lib/apt/lists/*

# Not shipped by dolfinx/dolfinx:nightly; everything else fem.py imports
# (dolfinx, ufl, mpi4py, petsc4py, numpy) is already on the base image.
RUN pip install --no-cache-dir gmsh pyyaml

WORKDIR /app
COPY . .

# CMD supplies default args; a bash script can override both, e.g.
# `docker run <image> other.msh weak`.
ENTRYPOINT ["python3", "fem.py"]
CMD ["mesh/rectangle-400x400.msh", "strong"]
