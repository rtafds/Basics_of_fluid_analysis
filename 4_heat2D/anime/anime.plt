if(exist("n")==0 || n<0) n=n0
time = sprintf("%.2f step",n)
set title time 
#filename = sprintf("Eta1_9_%05d_p0000.txt", n)
filename = sprintf("temp_%05d.dat", n)
splot filename using 1:2:3 with pm3d

n=n+dn
if( n<n1)reread
undefine n
