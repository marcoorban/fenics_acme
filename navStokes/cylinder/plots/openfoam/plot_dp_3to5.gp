set terminal pngcairo enhanced font "Helvetica,14" size 1400,1000
set border linewidth 1.2
set tics in
set mxtics 2
set mytics 2
set grid xtics ytics mxtics mytics lc rgb '#e0e0e0' lw 1, lc rgb '#eeeeee' lw 1

# Okabe-Ito colorblind-safe palette, same convention as plots.gp/plot2d-3.gp/plot_re100_3to5.gp.
set style line 6 lc rgb '#0072B2' lw 2.4 dt 1   # delta_p(t) -- blue
set style line 7 lc rgb '#D55E00' lw 3.0 dt 1   # highlighted period -- vermillion

set output "plots/dp_3to5.png"
set xlabel "{/Italic t}"
set ylabel "{/Symbol D}{/Italic p}"
set xrange [3:5]
set yrange [2.35:2.55]
set key off
set title "Pressure difference (front - rear stagnation), Re=100, t=3-5s"

# Half-period samples: delta_p = p(front) - p(rear), evaluated at t0 + 1/(2*f),
# f = f(Cl) (period T_cl = 1/f = 0.332523s from Cl's own zero-crossings, mean
# of 5 periods), t0 = each Cl local max in this window (parabolic-vertex
# refined). See cylinder/log.md's delta_p_at_half_period writeup for the
# method/rationale (DFG benchmark definition).
$dp_half_period << EOD
3.4847 2.459297
3.8166 2.464694
4.1490 2.466396
4.4817 2.468752
4.8143 2.469314
EOD

# Period boundaries and highlight below use f(Cl) -- same positive-going
# zero-crossings used for the delta_p half-period samples above (T_cl =
# 0.332523s, mean of 5 periods). Note: since delta_p/Cd complete two cycles
# per Cl period, one Cl-period-wide highlight will visibly span 2 delta_p
# oscillations -- that's expected here, not a bug (see the T_cl/T_cd~2 ratio
# noted in the plot_re100_3to5.gp / log.md writeups if a per-delta_p-cycle
# highlight is ever wanted instead). The last full Cl period in this window,
# [4.5664, 4.8992], is highlighted in vermillion over the base curve.
period_bounds = "3.2366 3.5687 3.9011 4.2336 4.5664 4.8992"
do for [x in period_bounds] {
    set arrow from x, graph 0 to x, graph 1 nohead lc rgb '#999999' dt 3 lw 1.2
}

plot "openFoam/postProcessing/pressureDropDict/3/p" using 1:($2-$3) with line ls 6 title "{/Symbol D}p(t)", \
     "< awk '!/^#/ && $1>=4.5664 && $1<=4.8992' openFoam/postProcessing/pressureDropDict/3/p" \
       using 1:($2-$3) with line ls 7 title "highlighted period", \
     2.4773 with lines dt 2 lw 1.5 lc rgb "black" title "{/Symbol D}p ref. (2.4773)", \
     $dp_half_period using 1:2 with points pt 7 ps 1.8 lc rgb "black" title "t_0+T/2 sample (t_0 = C_L local max)"
