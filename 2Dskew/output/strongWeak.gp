set terminal pngcairo size 800,600
set output "strongWeak.png"

set xrange [-0.05:1.05]
set ylabel "u^h"
set key outside bottom horizontal
set title "Weak and strong BCs, y=1.5"

plot 	 "solutionStrongWeak_rectangle-10x10.dat" w l lw 2 lc "blue"  title "10x10"  , \
	 "solutionStrongWeak_rectangle-20x20.dat" w l lw 2 lc "orange"  title "20x20"  , \
	 "solutionStrongWeak_rectangle-40x40.dat" w l lw 2 lc "green"  title "40x40"  , \
	 "solutionStrongWeak_rectangle-80x80.dat" w l lw 2 lc "violet"  title "80x80"  , \
	 "solutionStrongWeak_rectangle-100x100.dat" w l lw 2 lc "black"  title "100x100"  , \
   "solutionStrongWeak_rectangle-400x400.dat" w l lw 2 lc "red" title "400x400", \

set output
system "xdg-open strongWeak.png"
