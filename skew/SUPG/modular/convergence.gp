set terminal pngcairo size 800,600
set output "convergence.png"
set logscale x
set logscale y
set xlabel "h"
set ylabel "L_2 error"
set title "Convergence analysis"
set xrange [0.00009:0.012]
set key right bottom
plot "convergence.dat" u 2:3 w lp lt 3 title "strong-weak", \
	 "" u 2:4 w lp lt 10 title "strong", \
	 "" u 2:5 w lp lt 1 title "weak"
