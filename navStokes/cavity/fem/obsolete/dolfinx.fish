#!/usr/bin/env fish

set meshFile cavity-2D.msh

# Run from this directory
cd (dirname (status -f)); or exit 1

# Copy the mesh file here
cp ../mesh/$meshFile .

set containerName dolfinx
docker rm -f $containerName 2>/dev/null; or true
docker run -it --rm --name $containerName -v $PWD:/project -w /project --user (id -u):(id -g) dolfinx/dolfinx:nightly
