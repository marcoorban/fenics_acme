#!/usr/bin/env fish

# Run from this directory
cd (dirname (status -f)); or exit 1

source ~/miniconda3/etc/fish/conf.d/conda.fish
conda activate fenics-env; or exit 1

set meshFile cavity_5432cells.msh
set meshStem (basename $meshFile .msh)
set reynolds_numbers 0.5 1 10 50 100 200 400

for re in $reynolds_numbers
    echo "=== Re = $re ==="
    python3 fem_strong.py $meshFile $re; or echo "fem_strong.py failed at Re=$re"
    python3 postprocess.py results/strong/$meshStem/Re_$re/u_p0_000000.vtu; or echo "postprocess.py failed for strong, Re=$re"

    python3 fem_weak.py $meshFile $re; or echo "fem_weak.py failed at Re=$re"
    python3 postprocess.py results/weak/$meshStem/Re_$re/u_p0_000000.vtu; or echo "postprocess.py failed for weak, Re=$re"

    python3 fem_weak_noSUPG.py $meshFile $re; or echo "fem_weak_noSUPG.py failed at Re=$re"
    python3 postprocess.py results/weak_noSUPG/$meshStem/Re_$re/u_p0_000000.vtu; or echo "postprocess.py failed for weak_noSUPG, Re=$re"

    python3 centerline_profile.py \
        results/strong/$meshStem/Re_$re/u_p0_000000.vtu \
        results/weak/$meshStem/Re_$re/u_p0_000000.vtu \
        results/weak_noSUPG/$meshStem/Re_$re/u_p0_000000.vtu; \
        or echo "centerline_profile.py failed at Re=$re"
end
