set xrange [0:3]
set yrange [0:3]
set view 40,60


set pm3d
unset surface


set xlabel 'X'
set ylabel 'Y'
set zlabel 'Deplacement Max'

#set view map

splot 'coord_stat.dat' u 1:2:3 w l

pause 10







