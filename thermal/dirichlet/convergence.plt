set terminal cairolatex pdf size 4in, 4in standalone color font ",9"
set output 'convergence.tex'

set xlabel '{\small $h$}'
set ylabel '{\small $\|u - u^h\|_{L^2}$}'
set xtics font ",5"
set ytics font ",5"
set ytics offset -1, 0

set logscale x 
set logscale y
set format y '{\tiny $10^{%L}$}'
set format x '{\tiny $10^{%L}$}'

set key off

set grid xtics mxtics ytics mytics lc rgb "#E0E0E0" dt 2 back

set lmargin 11
set rmargin 2
set tmargin 2
set bmargin 4

### --- Slope reference triangle ---
# Pick an anchor point (x0, y0) near your data where you want the triangle to sit,
# and a horizontal "run" factor (in log space) to control its size.
x0 = 0.04          # left x-coordinate of triangle (adjust to sit near your data)
y0 = 1e-4          # bottom y-coordinate of triangle (adjust to sit near your data)
run = 1.6          # how many x-fold the triangle spans horizontally (size control)
slope = 2.0        # your convergence order

x1 = x0 * run
y1 = y0 * (run**slope)   # rise determined by slope, so hypotenuse matches slope exactly

# Triangle vertices: (x0,y0) -> (x1,y0) -> (x1,y1) -> back to (x0,y0)
set object 1 polygon from x0,y0 to x1,y0 to x1,y1 to x0,y0
set object 1 fc rgb "gray40" fillstyle empty border lc rgb "gray20" lw 1.2

# Slope label, placed just right of the vertical leg
set label 1 '{\tiny $2$}' at x1*1.15, (y0*y1)**0.5 textcolor rgb "gray20"
set label 2 '{\tiny $1$}' at 5.5e-2, 7.8e-5 textcolor rgb "gray20"

plot 'convergence.dat' using 1:2 with lp lc "dark-magenta" pointsize 0.8 linewidth 1.5

set output
