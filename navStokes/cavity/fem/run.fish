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

# Built once here and passed to the solver script, instead of having
# it rebuild the same path from mesh/Re/its own name -- the two can't
# drift out of sync if there's only one place that constructs it.
set meshLabel (basename $meshFile .msh)
set solverScriptStem (path basename -E $solverScript)
set resultsDir results/Re_$reynolds/$meshLabel/$solverScriptStem

python3 $solverScript $meshFile $reynolds $resultsDir; or begin
    echo "$solverScript failed"
    exit 1
end

python3 plot_fields.py $resultsDir; or echo "plot_fields.py failed"
python3 plot_centerlines.py $resultsDir $reynolds; or echo "plot_centerlines.py failed"
python3 plot_peclet.py $meshFile $resultsDir; or echo "plot_peclet.py failed"

echo "Done -- see $resultsDir/fields.png, $resultsDir/centerlines.png, and $resultsDir/peclet.png"
