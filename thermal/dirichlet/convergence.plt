set terminal svg size 600,400 dynamic enhanced font 'Arial,12';
set output 'convergence.svg';

set title "Solution convergence"
set key right top
set xlabel 'element size h'
set ylabel 'L2 error'
set logscale x 
set logscale y
set format y "10^{%L}"
set format x "10^{%L}"
set key off
set grid mxtics mytics 
set grid xtics mxtics ytics mytics lc rgb "#E0E0E0" dt 2 back

set label 1 "slope=2.00" at graph 0.42, 0.40 front textcolor "dark-magenta"
#set label 2 "Euler, m=-2.24" at graph 0.60, 0.78 front textcolor "dark-spring-green"

plot 'convergence.dat' skip 1 using 1:2 with lp lc "dark-magenta"
