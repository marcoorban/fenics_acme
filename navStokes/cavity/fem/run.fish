#!/usr/bin/env fish

# Run from this directory regardless of where it's invoked from.
cd (dirname (status -f)); or exit 1

source ~/miniconda3/etc/fish/conf.d/conda.fish
conda activate fenics-env; or exit 1

if test (count $argv) -lt 2
    echo "Usage: run.fish <mesh.msh> <reynolds> [solver_script]"
    exit 1
end

set meshFile $argv[1]
set reynolds $argv[2]
set solverScript fem_newton.py
if test (count $argv) -ge 3
    set solverScript $argv[3]
end

# Must match the solver script's own results_dir (mesh_label =
# Path(mesh).stem, Re folder = f"Re_{reynolds:g}") -- both sides just
# take the CLI value as-is, so pass reynolds the same way each time
# (e.g. "400", not "400.0") to keep the two in sync.
set meshLabel (basename $meshFile .msh)
set resultsDir results/Re_$reynolds/$meshLabel

python3 $solverScript $meshFile $reynolds; or begin
    echo "$solverScript failed"
    exit 1
end

python3 plot_fields.py $resultsDir; or echo "plot_fields.py failed"
python3 plot_centerlines.py $resultsDir $reynolds; or echo "plot_centerlines.py failed"

echo "Done -- see $resultsDir/fields.png and $resultsDir/centerlines.png"
