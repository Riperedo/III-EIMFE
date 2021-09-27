set title "Structure factor"
set xrange [0.5:0.7]
set yrange [1e-6:10]
#set key bottom right
set logscale y
set format y "10^{%L}"
set xlabel "{/Symbol f}"
plot 'DAT\arrest_diagram.dat' w l lw 3 t "{/Symbol r} = 0.1", # "DAT\\diagrama de arresto.dat"