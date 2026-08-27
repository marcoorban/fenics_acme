set terminal pngcairo size 800,600

set output "results_gamma_pos1.png"

set ylabel "u^h"
set xlabel "x"
set xrange [0.7:1.00]
set yrange [0:1.3]
set key left bottom

plot "output/8-line.dat" u 1:2 w l title "8 elements" lw 2 lc "blue", \
     "output/16-line.dat" u 1:2 w l title "16 elements" lw 2 lc "violet", \
     "output/32-line.dat" u 1:2 w l title "32 elements" lw 2 lc "red", \
     "output/64-line.dat" u 1:2 w l title "64 elements" lw 2 lc "forest-green", \
     "output/128-line.dat" u 1:2 w l title "128 elements" lw 2 lc "black", \

