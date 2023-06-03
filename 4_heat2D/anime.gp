reset 
set nokey
set size square
set pm3d map
set contour
set ticslevel 0
unset colorbox
set cntrparam levels incremental -1,0.01,0
set xr[0:20]
set yr[0:20]
set palette defined (-1 "white", 0 "white")
set cbrange[-1:0]

set term gif animate
set output "2D-cavity_unsteady-anime.gif"

n0=00000
n1=00400
dn=00005

load "anime.plt"
