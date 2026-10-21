# Plot u along the vertical centerline (x=0.5) and v along the horizontal
# centerline (y=0.5) for all three BC methods (strong, weak, weak_noSUPG),
# overlaid with the matching Ghia, Ghia & Shin (1982) benchmark columns
# (see ghia_re below) -- these are exactly the two profiles Ghia
# published: Table I is u along the VERTICAL centerline, Table II is v
# along the HORIZONTAL one (see centerline_profile.py, which samples and
# writes both lines).
#
# Ghia overlay: pick the Reynolds number with ghia_re (must be one Ghia
# actually tabulated: 100, 400, 1000, 3200, 5000, 7500, 10000). lid_speed
# is this solver's U_lid at that Re (U_lid = Re*nu/L, nu=0.01, L=1),
# needed to undo Ghia's lid-speed normalization (u*=u/U_lid) since this
# solver's raw output is dimensional, not normalized -- at Re=100 this is
# exactly 1.0, so the Ghia points need no scaling there. Using this
# script on a datafile from a different Re than ghia_re will still plot,
# but the Ghia overlay won't correspond to the same case.
#
# Usage (mesh_label defaults to $CAVITY_MESH, the same env var the
# solver scripts read -- set it once in your shell instead of retyping
# the path every time; re_label has no env-var equivalent since it
# varies per solver invocation, not per sweep script):
#   export CAVITY_MESH=normal
#   gnuplot -e "re_label='100'; ghia_re=100" plot_centerline.gp
# or override the path outright:
#   gnuplot -e "datafile='results/Re_100/normal/dat/centerline.dat'; ghia_re=100" plot_centerline.gp
#
# datafile (horizontal) columns (see centerline_profile.py):
#   1=x 2=u_strong 3=v_strong 4=div_strong 5=u_weak 6=v_weak 7=div_weak
#   8=u_noSUPG 9=v_noSUPG 10=div_noSUPG
# centerline_vertical.dat columns: 1=y, 2-10 as above but along x=0.5

if (!exists("mesh_label")) mesh_label = system("echo -n $CAVITY_MESH")
if (strlen(mesh_label) == 0) mesh_label = "normal"
if (!exists("re_label")) re_label = "100"
if (!exists("datafile")) datafile = "results/Re_" . re_label . "/" . mesh_label . "/dat/centerline.dat"
if (!exists("ghia_re")) ghia_re = 100
vertical_datafile = datafile[1:strlen(datafile)-4] . "_vertical.dat"
outfile = datafile[1:strlen(datafile)-4] . ".png"
lid_speed = ghia_re * 0.01  # U_lid = Re*nu/L, this solver's nu=0.01, L=1
if (ghia_re == 100)   { ghia_col = 2 }
if (ghia_re == 400)   { ghia_col = 3 }
if (ghia_re == 1000)  { ghia_col = 4 }
if (ghia_re == 3200)  { ghia_col = 5 }
if (ghia_re == 5000)  { ghia_col = 6 }
if (ghia_re == 7500)  { ghia_col = 7 }
if (ghia_re == 10000) { ghia_col = 8 }

set terminal pngcairo noenhanced size 1400,650 font ",11"
set output outfile

# Colorblind-safe palette (Okabe-Ito), one dash pattern per method so
# lines stay distinguishable even in grayscale printouts.
set style line 1 lc rgb '#0072B2' lw 2.5 dashtype 1   # strongBC: solid blue
set style line 2 lc rgb '#D55E00' lw 2.5 dashtype 2   # weakBC w/ SUPG (fem_weak_SUPG.py): dashed vermillion
set style line 3 lc rgb '#009E73' lw 2.5 dashtype '.-' # weakBC (fem_weak.py): dot-dash green
set style line 4 lc rgb '#000000' pt 6 ps 1.3          # Ghia et al. (1982): open circles

set multiplot layout 1,2 title datafile
set grid

set xlabel "y"
set ylabel "u (vertical centerline, x=0.5)"
plot vertical_datafile using 1:2 with lines ls 1 title "strongBC", \
     vertical_datafile using 1:5 with lines ls 2 title "weakBC w/ SUPG", \
     vertical_datafile using 1:8 with lines ls 3 title "weakBC", \
     "../benchmarks/ghia1982_table1_u_vertical.dat" using 1:(column(ghia_col)*lid_speed) with points ls 4 title sprintf("Ghia et al. (1982), Re=%d", ghia_re)

set xlabel "x"
set ylabel "v (horizontal centerline, y=0.5)"
plot datafile using 1:3 with lines ls 1 title "strongBC", \
     datafile using 1:6 with lines ls 2 title "weakBC w/ SUPG", \
     datafile using 1:9 with lines ls 3 title "weakBC", \
     "../benchmarks/ghia1982_table2_v_horizontal.dat" using 1:(column(ghia_col)*lid_speed) with points ls 4 title sprintf("Ghia et al. (1982), Re=%d", ghia_re)

unset multiplot
print "Wrote ".outfile
