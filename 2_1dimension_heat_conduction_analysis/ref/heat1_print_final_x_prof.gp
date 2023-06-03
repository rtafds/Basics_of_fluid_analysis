#
# heat1_print_final_x_prof.gp
#
# final temperature profile at y=0.5 as a function of x
#
set xrange [0:32]
set yrange [0:0.5]
set xlabel 'i'
set ylabel 'temp at y=0.5'
plot '../data/temp.final_profile_x' w lp 
pause -1