set style line 1 linecolor "red" linetype 1 linewidth 2 pointtype 1 pointsize 2
set style line 2 linecolor "magenta" linetype 2 linewidth 2 pointtype 2 pointsize 2
set style line 3 linecolor "blue" linetype 1 linewidth 2 pointtype 8 pointsize 2
set style line 4 linecolor "dark-green" linetype 1 linewidth 2 pointtype 6 pointsize 2
set style line 5 linecolor "gold" linetype 1 linewidth 2 pointtype 3 pointsize 2
set style line 6 linecolor "black" linetype 1 linewidth 2 pointtype 3 pointsize 2


set yrange[0.15:1.05]
plot "data_Qt25.txt" using 1:2 with lines ls 2 title "q_t = 25 m^3/d", "data_default.txt" using 1:2 with lines ls 1 title "q_t = 50 m^3/d", "data_Qt75.txt" using 1:2 with lines ls 3 title "q_t = 75 m^3/d",
set terminal pngcairo enhanced
set output "Resultado4.png"
set grid xtics ytics lw 2
set grid
set xlabel "x (m)"
set ylabel "Saturação da água"
set key right top
replot
set terminal wxt
set output

