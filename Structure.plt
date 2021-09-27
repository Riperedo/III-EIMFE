set title "Structure factor"
set yrange [0:3]
set xrange [0:20]
#set key bottom right

set xlabel "k"
plot 'DAT\S0p1.dat' w l lw 3 t "{/Symbol r} = 0.1", 'DAT\S0p2.dat' w l lw 3 t "{/Symbol r} = 0.2", 'DAT\S0p3.dat' w l lw 3 t "{/Symbol r} = 0.3", 'DAT\S0p4.dat' w l lw 3 t "{/Symbol r} = 0.4", 'DAT\S0p5.dat' w l lw 3 t "{/Symbol r} = 0.5", 'DAT\S0p6.dat' w l lw 3 t "{/Symbol r} = 0.6", 'DAT\S0p7.dat' w l lw 3 t "{/Symbol r} = 0.7", 'DAT\S0p8.dat' w l lw 3 t "{/Symbol r} = 0.8", 'DAT\S0p9.dat' w l lw 3 t "{/Symbol r} = 0.9"