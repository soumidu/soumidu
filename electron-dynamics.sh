set terminal pngcairo size 750,550 enhanced color font 'Times-Roman,25' lw 1.8
set output 'spm.png'
set key at 270,0.148
set key box lt 1 lc rgb "white" lw 0
set key width 1
set key right top
set key spacing 1.2

set ylabel '{/Times-Italic z} [au]' font 'Times-Roman,30'
set xlabel '{/Times-Italic t/T}' font 'Times-Roman,30'

set ytics (-280,0,280) font 'Arial,25'

set xtics ('0' 0,'0.5' 54.165,'1' 108.33, '1.5' 162.495 ,'2' 216.66,'2.5' 270.825) font 'Arial,25'


set object circle at 38.7,-3 size 0.5 fillcolor rgb "red" fillstyle solid

set label "{/Times-Italic z}_0" at 37.5,-45 font 'Times-Roman,27'
set label "{/Times-Italic r}_e" at 7.13,-45 font 'Times-Roman,27'
set label "{/Times-Italic r}_q" at 31.52,-45 font 'Times-Roman,27'
set label "(b)" at -250,0.125 font 'Times-Roman,27'

set title "{/Symbol e} = 0.5" font 'Times-Roman,29'

unset label 
unset key
unset title 
unset object
p [-40:300][-280:280] "traj_spm.txt" u 2:($1/-1) w l lt 2 lw 2 lc rgb "black" title "{/Times-Italic F_x(t)}","traj_spm_08.txt" u 2:($1/-1)  w l lt 2 lw 2 lc rgb "red" title "{/Times-Italic F_y(t)}","traj_spm_-08.txt" u 2:($1/-1) w l lt 2 lw 2 lc rgb "blue" title "{/Times-Italic F(t)}","laser.txt" u 2:($1*150) w l dt 2 lw 2 lc rgb "#005A32" title "{/Times-Italic F_x(t)}"

