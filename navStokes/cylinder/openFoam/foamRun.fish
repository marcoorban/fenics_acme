#!/usr/bin/env fish

# Run from this directory
cd (dirname (status -f)); or exit 1
bash -lc "./Allclean"
docker run --rm -v $PWD:/project -w /project --user (id -u):(id -g) -e HOME=/tmp microfluidica/openfoam:org bash -lc "set -o pipefail; decomposePar && mpirun -np 16 --use-hwthread-cpus foamRun -solver incompressibleFluid -parallel | tee log.foamRun && reconstructPar -latestTime && foamPostProcess -solver incompressibleFluid -func forceCoeffsDict -latestTime"
