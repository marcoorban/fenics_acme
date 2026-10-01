set terminal pngcairo enhanced font "Helvetica,14" size 1400,1000
set border linewidth 1.2
set tics in
set mxtics 2
set mytics 2
set grid xtics ytics mxtics mytics lc rgb '#e0e0e0' lw 1, lc rgb '#eeeeee' lw 1

# Okabe-Ito colorblind-safe palette. Re=100/400/1000/3200/5000 each get their
# own color, shared between a Ghia marker and its matching OpenFOAM line --
# color carries the Re, open-marker-vs-solid-line carries benchmark-vs-
# simulation. Dashed/dotted lines were dropped: with 5 overlapping curves,
# telling dash patterns apart at crossings was harder to read than color.
set style line 1 lc rgb '#0072B2' pt 4  ps 1.7 lw 1.4   # Re=100  -- open square, blue
set style line 2 lc rgb '#D55E00' pt 6  ps 1.7 lw 1.4   # Re=400  -- open circle, vermillion
set style line 3 lc rgb '#009E73' pt 8  ps 1.7 lw 1.4   # Re=1000 -- open triangle up, bluish green
set style line 4 lc rgb '#CC79A7' pt 10 ps 1.7 lw 1.4   # Re=3200 -- open triangle down, reddish purple
set style line 5 lc rgb '#E69F00' pt 12 ps 1.7 lw 1.4   # Re=5000 -- open diamond, orange

set style line 6  lc rgb '#0072B2' lw 2.4 dt 1   # OpenFOAM Re=100  -- blue
set style line 7  lc rgb '#D55E00' lw 2.4 dt 1   # OpenFOAM Re=400  -- vermillion
set style line 8  lc rgb '#009E73' lw 2.4 dt 1   # OpenFOAM Re=1000 -- bluish green
set style line 9  lc rgb '#CC79A7' lw 2.4 dt 1   # OpenFOAM Re=3200 -- reddish purple
set style line 10 lc rgb '#E69F00' lw 2.4 dt 1   # OpenFOAM Re=5000 -- orange

set output "plots/forceCoeffs.png"
set xlabel "{/Italic t}"
set ylabel "{/Italic Cf}"
set xrange [0:8]
set yrange [-0.8:3.2]
set key off
#set label 1 "Re = 100" at graph 0.05, 0.90

plot "openFoam/postProcessing/forceCoeffsDict/0/forceCoeffs.dat" using 1:3 with line ls 6 title "Cd", \
     "openFoam/postProcessing/forceCoeffsDict/0/forceCoeffs.dat" using 1:4 with line ls 7 title "Cl", \
     2.9483 with lines dt 2 lw 1.5 lc rgb "black" title "Cd, 9a", \
     0.4651 with lines dt 2 lw 1.5 lc rgb "black" title "Cl, 9a"
