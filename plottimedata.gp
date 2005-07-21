set xlabel "n";
set ylabel "m";
set zlabel "sec/iter";
##--
set dgrid3d 5,5,3;
set pm3d;

splot "timedata.dat" using 1:2:(1/$3) with lp title "las2";
