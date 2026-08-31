set terminal pngcairo size 800, 600 enhanced
set output "convergence.png"

set logscale x
set logscale y 
set xlabel "h" font "IBM Plex Sans Italic"
set ylabel "|| u - u^h ||_{L^2}" font "IBM Plex Sans Italic"
set format y "10^{%T}"; set tics (0.1, 0.01, 0.001, 0.001) font "IBM Plex Sans,8"
set format x "10^{%T}"; set tics (0.1, 0.01, 0.001, 0.001) font "IBM Plex Sans,8"
set xrange [0.0005:0.6]
set grid 
set grid mytics
set grid mxtics

plot "convergence.dat" u 2:3 w lp ls 1 title "{/Symbol g}=-1", \
      ""               u 2:4 w lp ls 2 title "{/Symbol g}=+1"


