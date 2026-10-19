#!/usr/bin/env fish

set -l script_dir (dirname (status --current-filename))
set -l mpirun /home/morban/miniconda3/envs/fenics-env/bin/mpirun
set -l python /home/morban/miniconda3/envs/fenics-env/bin/python
set -l num_procs 16
set -l mesh_file quad2D_2504cells.msh
set -l mesh_name (string replace -r '\.msh$' '' $mesh_file)

set -l results_dir $script_dir/../../results
set -l comparison_csv $results_dir/comparison.csv
mkdir -p $results_dir
echo "solver,weakStrong,case,mesh,mesh_elements,mesh_nodes,num_threads,wall_time,Umax,Pmax,CD_max,CL_max" >$comparison_csv

for bc in strong weak
    for case in 1 2 3
        echo "=== Running case $case ($bc) ==="
        set -l start (date +%s)
        set -l output ($mpirun -n $num_procs $python $script_dir/main.py --case $case --bc $bc --mesh $mesh_file)
        set -l exit_code $status
        set -l wall_time (math (date +%s) - $start)

        if test $exit_code -ne 0
            echo "!!! case $case ($bc) failed with exit code $exit_code -- skipping its comparison row !!!"
            continue
        end

        # main.py's last stdout line: num_cells,num_nodes,Umax,Pmax,CD_max,CL_max
        set -l summary (string split , -- $output[-1])
        set -l mesh_elements $summary[1]
        set -l mesh_nodes $summary[2]
        set -l u_max $summary[3]
        set -l p_max $summary[4]
        set -l cd_max $summary[5]
        set -l cl_max $summary[6]

        echo "chorin,$bc,$case,$mesh_name,$mesh_elements,$mesh_nodes,$num_procs,$wall_time,$u_max,$p_max,$cd_max,$cl_max" >>$comparison_csv
    end
end

echo "All 6 runs finished. Comparison table: $comparison_csv"
