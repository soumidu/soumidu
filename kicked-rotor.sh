set terminal pngcairo size 750,750 enhanced color font 'Times-Roman,27' lw 1.8
set output 'kicked-rotor.png'
set key at 12e14,30
set key box lt 1 lc rgb "white" lw 0
set key width 1
set key right top
set key spacing 1.2

set xlabel 'q' font 'Times-Roman,28'
set ylabel 'p' font 'Times-Roman,28'

set yrange [-3.14:3.14]
set xrange [0:6]
set ytics (-3,0,3) font 'Arial,30'
set xtics (0,3,6) font 'Arial,30'

 p "sm.dat" u 1:2:-1 lc var ps 0.5 pt 7


