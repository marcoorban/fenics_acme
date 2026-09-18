set terminal pngcairo enhanced font "Helvetica,14" size 1400,1000
set border linewidth 1.2
set tics in
set mxtics 2
set grid xtics mxtics lc rgb '#e0e0e0' lw 1

# Okabe-Ito colorblind-safe palette, same convention as the other cylinder plots.
# Cl is drawn uniformly faded throughout (alpha channel on the pngcairo-
# supported '#aarrggbb' color form, no exceptions/highlights) since it's now
# just visual context for delta_p, not a primary series -- and its own axis
# is hidden below. Two solid-color point markers (see $cl_max_pointers below)
# pin the specific Cl-max values of interest instead of a full tics/label axis.
set style line 6 lc rgb '#0072B2' lw 2.4 dt 1     # delta_p(t) -- blue, right axis
set style line 7 lc rgb '#C0D55E00' lw 2.4 dt 1   # Cl(t) -- faded vermillion, left axis (hidden)

set output "plots/cl_dp_3to5.png"
set xlabel "{/Italic t}"
set xrange [3:5]
set key off
set title "C_L and {/Symbol D}p superimposed, Re=100, t=3-5s"

# Left axis (Cl) hidden entirely -- no label, no tics -- since only its shape/
# timing matters here, not its values (delta_p's axis on the right is what
# the reader actually reads numbers off).
unset ylabel
unset ytics
set yrange [-1.1:1.1]

set y2label "{/Symbol D}{/Italic p}" tc rgb '#0072B2'
set y2range [2.35:2.55]
set y2tics nomirror tc rgb '#0072B2'

# Two Cl-max pointers (genuine local maxima of Cl, parabolic-vertex refined)
# bounding one full Cl period: t=4.3154 and t=4.6480 (spacing 0.3326s ~ T_cl).
# No Cl-min pointer per request.
$cl_max_pointers << EOD
4.3154 0.9538957
4.6480 0.9527795
EOD

# The single delta_p value between them: p(front)-p(rear) at the exact
# midpoint of the two Cl-max pointers above (t=4.4817, smack in the middle
# of the two peaks) -- see cylinder/log.md's delta_p_at_half_period writeup
# for the method/rationale.
$dp_between << EOD
4.4817 2.468795
EOD

# Period boundaries: Cl's own peak times (NOT zero-crossings) -- each period
# starts and ends at the same y-value (a Cl maximum), rather than at Cl=0.
# Kept as faint markers only (no colored highlight -- Cl stays uniformly
# transparent everywhere per request).
period_bounds = "3.3184 3.6504 3.9828 4.3154 4.6480 4.9809"
do for [x in period_bounds] {
    set arrow from x, graph 0 to x, graph 1 nohead lc rgb '#999999' dt 3 lw 1.2
}

plot "openFoam/postProcessing/forceCoeffsDict/3/forceCoeffs.dat" using 1:4 axes x1y1 \
       with line ls 7 title "Cl", \
     "openFoam/postProcessing/pressureDropDict/3/p" using 1:($2-$3) axes x1y2 \
       with line ls 6 title "{/Symbol D}p", \
     2.4773 axes x1y2 with lines dt 2 lw 1.5 lc rgb "black" notitle, \
     $cl_max_pointers using 1:2 axes x1y1 with points pt 7 ps 1.8 lc rgb '#D55E00' notitle, \
     $dp_between using 1:2 axes x1y2 with points pt 7 ps 1.8 lc rgb "black" notitle
