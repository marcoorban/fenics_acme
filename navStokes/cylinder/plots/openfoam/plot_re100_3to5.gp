set terminal pngcairo enhanced font "Helvetica,14" size 1400,1000
set border linewidth 1.2
set tics in
set mxtics 2
set mytics 2
set grid xtics ytics mxtics mytics lc rgb '#e0e0e0' lw 1, lc rgb '#eeeeee' lw 1

# Okabe-Ito colorblind-safe palette, same convention as plots.gp/plot2d-3.gp.
set style line 6 lc rgb '#0072B2' lw 2.4 dt 1   # Cd -- blue
set style line 7 lc rgb '#D55E00' lw 2.4 dt 1   # Cl -- vermillion

set output "plots/forceCoeffs_3to5.png"
set xlabel "{/Italic t}"
set ylabel "{/Italic C_D}, {/Italic C_L}"
set xrange [3:5]
#set key top left box opaque
set key off
set title "Drag and lift coefficients, Re=100, t=3-5s"

plot "openFoam/postProcessing/forceCoeffsDict/3/forceCoeffs.dat" using 1:3 with line ls 6 title "Cd", \
     "openFoam/postProcessing/forceCoeffsDict/3/forceCoeffs.dat" using 1:4 with line ls 7 title "Cl", \
     3.2232 with lines dt 2 lw 1.5 lc rgb "black" title "Cd ref. (3.2232)", \
     0.9830 with lines dt 2 lw 1.5 lc rgb "black" title "Cl ref. (0.9830)"
