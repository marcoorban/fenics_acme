# Lid-driven cavity, Re=100 -- validation against Ghia, Ghia & Shin (1982).
# Ghia's tabulated data (17 discrete grid points) is drawn as open markers;
# the OpenFOAM profile (200-point sampleDict line) is drawn as a solid line
# -- reference vs. present-result is the standard convention for this kind
# of validation figure.

set terminal pngcairo enhanced font "Helvetica,14" size 1200,900
set border linewidth 1.2
set tics in
set mxtics 2
set mytics 2
set grid xtics ytics mxtics mytics lc rgb '#e0e0e0' lw 1, lc rgb '#eeeeee' lw 1

# Okabe-Ito colorblind-safe palette: black (reference) + blue (present result)
set style line 1 lc rgb '#000000' pt 6 ps 1.3 lw 1.4   # Ghia (1982) -- markers only
set style line 2 lc rgb '#0072B2' lt 1 lw 2.4           # OpenFOAM   -- line only

set output "plots/u.png"
set xlabel "{/Italic y}"
set ylabel "{/Italic u}"
set xrange [0:1]
set key top left box opaque
#set label 1 "Re = 100" at graph 0.05, 0.90

plot "benchmarks/ghia1982_table1_u_vertical.dat" using 1:2 with points ls 1 title "Ghia et al. (1982)", \
     "benchmarks/ghia1982_table1_u_vertical.dat" using 1:3 with points ls 1 title "", \
     "benchmarks/ghia1982_table1_u_vertical.dat" using 1:4 with points ls 1 title "", \
     "openFoam-results/re100/verticalCenterline.xy" using 1:($2/1) with lines ls 2 title "OpenFOAM", \
     "openFoam-results/re400/verticalCenterline.xy" using 1:($2/4) with lines ls 2 title "", \
     "openFoam-results/re1000/verticalCenterline.xy" using 1:($2/10) with lines ls 2 title "", \

set output "plots/v.png"
set xlabel "{/Italic x}"
set ylabel "{/Italic v}"
set xrange [0:1]
set key top right box opaque
#set label 1 "Re = 100" at graph 0.05, 0.90

plot "benchmarks/ghia1982_table2_v_horizontal.dat" using 1:2 with points ls 1 title "Ghia et al. (1982)", \
     "benchmarks/ghia1982_table2_v_horizontal.dat" using 1:3 with points ls 1 title "", \
     "benchmarks/ghia1982_table2_v_horizontal.dat" using 1:4 with points ls 1 title "", \
     "openFoam-results/re100/horizontalCenterline.xy" using 1:($3/1) with lines ls 2 title "OpenFOAM", \
     "openFoam-results/re400/horizontalCenterline.xy" using 1:($3/4) with lines ls 2 title "", \
     "openFoam-results/re1000/horizontalCenterline.xy" using 1:($3/10) with lines ls 2 title "", \
