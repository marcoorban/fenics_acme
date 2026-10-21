# Plot u, v and div(u) along the cavity's horizontal centerline for all
# three BC methods (strong, weak, weak_noSUPG), plus u along the vertical
# centerline, stacked into one image. The v (horizontal) and u (vertical)
# panels are overlaid with the matching Ghia, Ghia & Shin (1982) Re=400
# benchmark columns -- Ghia's Table I tabulates u along the VERTICAL
# centerline (x=0.5) and Table II tabulates v along the HORIZONTAL one
# (y=0.5), not both along either, so each panel needs its own line's data
# (see centerline_profile.py, which now samples and writes both).
#
# Ghia-overlay scope: Re=400 only, for now -- lid_speed below is this
# solver's U_lid at Re=400 (U_lid = Re*nu/L = 400*0.01/1 = 4), needed to
# undo Ghia's lid-speed normalization (u*=u/U_lid) since this solver's
# raw output is dimensional, not normalized. Using this script on a
# datafile from a different Re will still plot, but the Ghia overlay
# will not correspond to the same case.
#
# Usage:
#   gnuplot -e "datafile='results/dat/cavity_5432cells_Re_400.dat'" plot_centerline.gp
#
# datafile (horizontal) columns (see centerline_profile.py):
#   1=x 2=u_strong 3=v_strong 4=div_strong 5=u_weak 6=v_weak 7=div_weak
#   8=u_noSUPG 9=v_noSUPG 10=div_noSUPG
# <same base name>_vertical.dat columns: 1=y, 2-10 as above but along x=0.5

if (!exists("datafile")) datafile = "results/dat/cavity_5432cells_Re_400.dat"
vertical_datafile = datafile[1:strlen(datafile)-4] . "_vertical.dat"
outfile = datafile[1:strlen(datafile)-4] . ".png"
lid_speed = 4.0  # U_lid at Re=400 -- see header note

set terminal pngcairo noenhanced size 1400,1200 font ",11"
set output outfile

# Colorblind-safe palette (Okabe-Ito), one dash pattern per method so
# lines stay distinguishable even in grayscale printouts.
set style line 1 lc rgb '#0072B2' lw 2.5 dashtype 1   # strong: solid blue
set style line 2 lc rgb '#D55E00' lw 2.5 dashtype 2   # weak: dashed vermillion
set style line 3 lc rgb '#009E73' lw 2.5 dashtype '.-' # weak_noSUPG: dot-dash green
set style line 4 lc rgb '#000000' pt 6 ps 1.3          # Ghia et al. (1982): open circles

set multiplot layout 2,2 title datafile
set grid

set xlabel "x"
set ylabel "u"
plot datafile using 1:2 with lines ls 1 title "strong", \
     datafile using 1:5 with lines ls 2 title "weak", \
     datafile using 1:8 with lines ls 3 title "weak_noSUPG"

set xlabel "x"
set ylabel "v"
plot datafile using 1:3 with lines ls 1 title "strong", \
     datafile using 1:6 with lines ls 2 title "weak", \
     datafile using 1:9 with lines ls 3 title "weak_noSUPG", \
     "../benchmarks/ghia1982_table2_v_horizontal.dat" using 1:($3*lid_speed) with points ls 4 title "Ghia et al. (1982), Re=400"

set xlabel "x"
set ylabel "div(u)"
plot datafile using 1:4 with lines ls 1 title "strong", \
     datafile using 1:7 with lines ls 2 title "weak", \
     datafile using 1:10 with lines ls 3 title "weak_noSUPG"

set xlabel "y"
set ylabel "u (vertical centerline, x=0.5)"
plot vertical_datafile using 1:2 with lines ls 1 title "strong", \
     vertical_datafile using 1:5 with lines ls 2 title "weak", \
     vertical_datafile using 1:8 with lines ls 3 title "weak_noSUPG", \
     "../benchmarks/ghia1982_table1_u_vertical.dat" using 1:($3*lid_speed) with points ls 4 title "Ghia et al. (1982), Re=400"

unset multiplot
print "Wrote ".outfile
