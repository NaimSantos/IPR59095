#set yrange[0.4:1.6]
plot "output_data.txt" using 1:2 notitle w l lw 2 lc "blue"
set terminal pngcairo enhanced
set output "pressao_no_espaco.png"
set grid xtics ytics lw 2
set grid
set xlabel "x (m)"
set ylabel "Pressao" rotate by 0
replot
set terminal wxt
set output


# with lines = w l
# linewidth = lw
# linecolor = lc
# dashtype = dt