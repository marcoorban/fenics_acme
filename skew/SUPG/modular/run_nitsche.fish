#!/usr/bin/env fish
# Run the "nitsche" solver image, forwarding mesh/solver args to fem.py
# inside the container. Nothing here requires an image rebuild:
#
#   - ./mesh is bind-mounted read-only over /app/mesh, so dropping a new or
#     edited .msh file in ./mesh makes it visible to the container
#     immediately, instead of only whatever was COPY'd in at build time;
#   - ./output is bind-mounted read-write over /app/output, so the .xdmf
#     results land back on the host instead of staying trapped inside the
#     container (which is removed on exit via --rm).
#
# Usage:
#   ./run_nitsche.fish [mesh_file] [solver] [extra fem.py args...]
#     mesh_file  path as seen inside the container (default: mesh/rectangle.msh)
#     solver     strong | weak | weakStrong    (default: strong)
#     extra args forwarded as-is to fem.py, e.g. --sample-line

set -l script_dir (dirname (status -f))

set -l mesh_file $argv[1]
if test -z "$mesh_file"
    set mesh_file "mesh/rectangle.msh"
end

set -l solver $argv[2]
if test -z "$solver"
    set solver "strong"
end

set -l extra_args $argv[3..-1]

mkdir -p $script_dir/output

docker run --rm \
    -v $script_dir/mesh:/app/mesh:ro \
    -v $script_dir/output:/app/output \
    nitsche $mesh_file $solver $extra_args
