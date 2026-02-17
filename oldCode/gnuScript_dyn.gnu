#set term gif animate
#set output 'gif_output.gif'

#view#
set view 40,60
#set view map

#label#
set xlabel 'X'
set ylabel 'Y'
set zlabel 'u(x,y)'

#range#
set xrange [0:3]
set yrange [0:3]
set zrange[0:20]
set cbrange[0:20]

#3d parameter#
set pm3d
unset surface

#constants#
filedata = 'coord_dyn.dat'
duree = 0
i = 0