#!/usr/bin/env fish

# Run from this directory
cd (dirname (status -f)); or exit 1

cp ../mesh/cavity.msh .

set meshFile "cavity.msh"

docker run --rm -v $PWD:/project -w /project --user (id -u):(id -g) -e HOME=/tmp microfluidica/openfoam:org bash -lc "gmshToFoam $meshFile &&
  foamDictionary -set 'entry0/empty/type=empty, entry0/wall/type=wall' \
    constant/polyMesh/boundary"
