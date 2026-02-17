c = 3 + i

set title sprintf('m = %d Time (s)= %.3f', m, duree)

splot "coord_dyn.dat" u 1:2:c w pm3d
pause 0.2

i = i + 1
duree = i*pas

if( i < imax ) reread
