#
# 失敗作
#

set xlabel 'j'              # x-axis
set ylabel 'temperature'    # y-axis
set xrange[0:21]     # j-grid min & max
set yrange[0.0:0.5]         # temperature min & max

pause 5
do for[i = 1: 50]{
    input = sprintf("./anime/temp_%05d.dat", i)
    set term postscript enhanced eps color
    plot input w lp
    pause 1
}

