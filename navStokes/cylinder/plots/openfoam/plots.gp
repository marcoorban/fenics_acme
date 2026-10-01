# DFG 2D-2 benchmark (Re=100, unsteady vortex shedding) validation plots.
# Cd/Cl time histories over the periodic tail (t>=3s, ~5 shedding cycles),
# following the project plotting convention (see navStokes/log.md #plotting):
# Okabe-Ito colorblind-safe palette, pngcairo, italicized physical-quantity
# axis labels. No tabulated benchmark *time series* exists to overlay as
# markers here (unlike cavity's Ghia profiles) -- DFG 2D-2 publishes only
# scalar reference ranges (St, c_D,mean, c_L,max), called out as horizontal
# reference lines instead.

set terminal pngcairo enhanced font "Helvetica,14" size 1000,600
set datafile separator whitespace

set xlabel "{/Italic t} [s]"
set grid

# --- Cd(t) ---
set output "plots/cd_vs_t.png"
set ylabel "{/Italic C_D}"
set xrange [3:5]
set key top right
set title "Drag coefficient, DFG 2D-2 (Re=100)"
plot "openFoam/postProcessing/forceCoeffsDict/0/forceCoeffs.dat" using 1:3 \
       with lines lw 2 lc rgb "#0072B2" title "simulation", \
     3.2 with lines dt 2 lw 1.5 lc rgb "black" title "benchmark ref. (mean {/Symbol \273} 3.2)"

# --- Cl(t) ---
set output "plots/cl_vs_t.png"
set ylabel "{/Italic C_L}"
set xrange [3:5]
set key top right
set title "Lift coefficient, DFG 2D-2 (Re=100)"
plot "openFoam/postProcessing/forceCoeffsDict/0/forceCoeffs.dat" using 1:4 \
       with lines lw 2 lc rgb "#0072B2" title "simulation", \
      1.0 with lines dt 2 lw 1.5 lc rgb "black" title "benchmark ref. (max {/Symbol \273} 1.0)", \
     -1.0 with lines dt 2 lw 1.5 lc rgb "black" notitle

# --- delta_p(t) ---
# The benchmark's delta_p is NOT a window statistic of this signal -- it's
# defined at one specific instant per cycle, t0+T/2 where t0 is a Cl local
# max (see plot_results.py's delta_p_at_half_period / log.md 2026-09-14
# correction entry). The two points below are exactly the samples that
# function computed; marked here rather than implying the whole oscillating
# trace itself is the benchmark-comparable quantity.
set output "plots/dp_vs_t.png"
set ylabel "{/Symbol D}{/Italic p}"
set xrange [3:5]
set key top right
set title "Pressure drop (front - rear stagnation), DFG 2D-2 (Re=100)"
$dp_half_period << EOD
4.481581 2.468752
4.814445 2.469314
EOD
plot "< awk '!/^#/{print $1, $2-$3}' openFoam/postProcessing/pressureDropDict/0/p" \
       using 1:2 with lines lw 2 lc rgb "#0072B2" title "{/Symbol D}p(t)", \
     $dp_half_period using 1:2 with points pt 7 ps 1.5 lc rgb "black" \
       title "benchmark-defined sample (t_0+T/2, t_0 = C_L local max)"
