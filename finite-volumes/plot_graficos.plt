plot "data.txt" using 1:2 title "Exato" w l lw 2 lc "forest-green", "data.txt" using 1:3 title "Numérico" w l lw 2 lc "blue" dt '-'
set terminal pngcairo
set output "Resultado.png"
set grid
#set title "Comparação das soluções"
set xlabel "x"
set ylabel "Phi(x)"
replot
set terminal wxt
set output


# with lines = w l
# linewidth = lw
# linecolor = lc
# dashtype = dt