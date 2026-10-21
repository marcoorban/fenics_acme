#!/usr/bin/env fish

# Run from this directory
cd (dirname (status -f)); or exit 1

source ~/miniconda3/etc/fish/conf.d/conda.fish
conda activate fenics-env; or exit 1

set meshFile cavity_1968cells.msh
set -gx CAVITY_MESH coarser  # read by fem_*.py, centerline_profile.py, plot_centerline.gp
set reynolds_numbers 0.5 1 10 50 100 200 400

for re in $reynolds_numbers
    echo "=== Re = $re ==="
    set reDir results/Re_$re/$CAVITY_MESH

    python3 fem_strong.py $meshFile $re; or echo "fem_strong.py failed at Re=$re"
    python3 postprocess.py $reDir/strong/u_p0_000000.vtu; or echo "postprocess.py failed for strong, Re=$re"

    python3 fem_weak_SUPG.py $meshFile $re; or echo "fem_weak_SUPG.py failed at Re=$re"
    python3 postprocess.py $reDir/weak_SUPG/u_p0_000000.vtu; or echo "postprocess.py failed for weak_SUPG, Re=$re"

    python3 fem_weak.py $meshFile $re; or echo "fem_weak.py failed at Re=$re"
    python3 postprocess.py $reDir/weak/u_p0_000000.vtu; or echo "postprocess.py failed for weak, Re=$re"

    python3 centerline_profile.py \
        $reDir/strong/u_p0_000000.vtu \
        $reDir/weak_SUPG/u_p0_000000.vtu \
        $reDir/weak/u_p0_000000.vtu; \
        or echo "centerline_profile.py failed at Re=$re"
end
