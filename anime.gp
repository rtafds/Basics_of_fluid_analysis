reset 
set nokey
set size square
set pm3d map
set contour
set ticslevel 0
unset colorbox
set cntrparam levels incremental -1,0.01,0
set xr[1:64]
set yr[1:64]
set palette defined (-1 "white", 0 "white")
set cbrange[-1:0]

set term gif animate
set output "2D-cavity_unsteady-anime.gif"

n0=100
n1=29900
dn=100

load "anime.plt"
