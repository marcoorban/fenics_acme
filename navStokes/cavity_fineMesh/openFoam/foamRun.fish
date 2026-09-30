#!/usr/bin/env fish

# Run from this directory
cd (dirname (status -f)); or exit 1

# Stop and remove any still-running instance of this script before starting
# a new one. Without this, re-running the script while a previous run is
# still going starts a SECOND container that tees into the same
# log.foamRun concurrently, interleaving/corrupting it -- caught 2026-09-18
# with three concurrent runs left going at once, all writing the same file.
set containerName cavity-foamrun
docker rm -f $containerName 2>/dev/null; or true

bash -lc "./Allclean"

# Clear the log file. tee would truncate it anyway at the start of a clean
# run, but this also wipes any interleaved leftovers from a prior run that
# just got force-stopped above before it could exit/close the file cleanly.
: > log.foamRun

docker run --rm --name $containerName -v $PWD:/project -w /project --user (id -u):(id -g) -e HOME=/tmp microfluidica/openfoam:org bash -lc "decomposePar && mpirun -np 8 foamRun -solver incompressibleFluid -parallel | tee log.foamRun && reconstructPar -latestTime && foamPostProcess -func sampleDict -latestTime"
