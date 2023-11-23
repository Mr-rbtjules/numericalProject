set xrange [0:6]
set yrange [0:10]
set view 40,60

set pm3d
unset surface

set xlabel 'X'
set ylabel 'Y'
set zlabel 'Deplacement Max'


set size 0.7, 1.4  
set origin 0.15, -0.2 

splot 'coord_stat.dat' u 1:2:3 w l

pause 10
