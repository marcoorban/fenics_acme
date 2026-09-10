#!/usr/bin/env fish

# Run from this directory
cd (dirname (status -f)); or exit 1
bash -lc "./Allclean"
docker run --rm -v $PWD:/project -w /project --user (id -u):(id -g) -e HOME=/tmp microfluidica/openfoam:org bash -lc "decomposePar && mpirun -np 8 foamRun -solver incompressibleFluid -parallel | tee log.foamRun && reconstructPar -latestTime && foamPostProcess -func sampleDict -latestTime"
