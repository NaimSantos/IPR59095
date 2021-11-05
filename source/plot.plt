plot "output_data_1.txt" using 1:2 notitle w l lw 2 lc "blue"
set terminal pngcairo enhanced
set output "pressao_por_spaco.png"
set grid xtics ytics lw 2
set grid
set xlabel "x (m)"
set ylabel "Pressao (kPa)"# rotate by 0
replot
set terminal wxt
set output


# with lines = w l
# linewidth = lw
# linecolor = lc
# dashtype = dt