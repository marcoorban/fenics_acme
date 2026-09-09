# Lid-driven cavity, Re=100 -- validation against Ghia, Ghia & Shin (1982).
# Ghia's tabulated data (17 discrete grid points) is drawn as open markers;
# the OpenFOAM profile (200-point sampleDict line) is drawn as a solid line
# -- reference vs. present-result is the standard convention for this kind
# of validation figure.

set terminal pngcairo enhanced font "Helvetica,14" size 1000,1000
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

set output "plots/u.png"
set xlabel "{/Italic y}"
set ylabel "{/Italic u}"
set xrange [0:1]
set key top left box opaque
#set label 1 "Re = 100" at graph 0.05, 0.90

plot "benchmarks/ghia1982_table1_u_vertical.dat" using 1:2 with points ls 1 title "Ghia et al. (1982)", \
     "benchmarks/ghia1982_table1_u_vertical.dat" using 1:3 with points ls 2 title "", \
     "benchmarks/ghia1982_table1_u_vertical.dat" using 1:4 with points ls 3 title "", \
     "benchmarks/ghia1982_table1_u_vertical.dat" using 1:5 with points ls 4 title "", \
     "benchmarks/ghia1982_table1_u_vertical.dat" using 1:6 with points ls 5 title "", \
     "openFoam-results/re100/verticalCenterline.xy" using 1:($2/1) with lines ls 6 title "OpenFOAM Re=100", \
     "openFoam-results/re400/verticalCenterline.xy" using 1:($2/4) with lines ls 7 title "OpenFOAM Re=400", \
     "openFoam-results/re1000/verticalCenterline.xy" using 1:($2/10) with lines ls 8 title "OpenFOAM Re=1000", \
     "openFoam-results/re3200/verticalCenterline.xy" using 1:($2/32) with lines ls 9 title "OpenFOAM Re=3200", \
     "openFoam-results/re5000/verticalCenterline.xy" using 1:($2/50) with lines ls 10 title "OpenFOAM Re=5000", \

set output "plots/v.png"
set xlabel "{/Italic x}"
set ylabel "{/Italic v}"
set xrange [0:1]
set key top right box opaque
#set label 1 "Re = 100" at graph 0.05, 0.90

plot "benchmarks/ghia1982_table2_v_horizontal.dat" using 1:2 with points ls 1 title "Ghia et al. (1982)", \
     "benchmarks/ghia1982_table2_v_horizontal.dat" using 1:3 with points ls 2 title "", \
     "benchmarks/ghia1982_table2_v_horizontal.dat" using 1:4 with points ls 3 title "", \
     "benchmarks/ghia1982_table2_v_horizontal.dat" using 1:5 with points ls 4 title "", \
     "benchmarks/ghia1982_table2_v_horizontal.dat" using 1:6 with points ls 5 title "", \
     "openFoam-results/re100/horizontalCenterline.xy" using 1:($3/1) with lines ls 6 title "OpenFOAM Re=100", \
     "openFoam-results/re400/horizontalCenterline.xy" using 1:($3/4) with lines ls 7 title "OpenFOAM Re=400", \
     "openFoam-results/re1000/horizontalCenterline.xy" using 1:($3/10) with lines ls 8 title "OpenFOAM Re=1000", \
     "openFoam-results/re3200/horizontalCenterline.xy" using 1:($3/32) with lines ls 9 title "OpenFOAM Re=3200", \
     "openFoam-results/re5000/horizontalCenterline.xy" using 1:($3/50) with lines ls 10 title "OpenFOAM Re=5000", \
