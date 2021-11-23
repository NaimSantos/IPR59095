plot "time_evolution_data.txt" using 1:2 title "x = 0 (m)" w l lw 2 lc "blue","time_evolution_data.txt" using 1:3 title "x = 0.1L (m)" w l lw 2 lc "red", "time_evolution_data.txt" using 1:4 title "x = 0.2L (m)" w l lw 2 lc "green", "time_evolution_data.txt" using 1:5 title "x = 0.9L (m)" w l lw 2 lc "orange"
set terminal pngcairo enhanced
set output "pressao_no_tempo.png"
set grid xtics ytics lw 2
set grid
set key bottom left
set xlabel "t (dias)"
set ylabel "Pressao (kPa)"# rotate by 0
replot
set terminal wxt
set output


# with lines = w l
# linewidth = lw
# linecolor = lc
# dashtype = dt