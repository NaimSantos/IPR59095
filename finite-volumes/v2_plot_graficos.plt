set encoding utf8
s(x) = (x>=0.6 && x<=0.8) ? 1 : 0
f(x) = exp(-200*(x - 0.3)**2) + s(x)
set xrange[0:1]
plot "data.txt" using 1:3 title "Numérico" w l lw 3 lc "blue" dt '_', f(x) title "Exato" w l lw 2 lc "red"
set terminal pngcairo enhanced
set output "Resultado.png"
set grid
#set title "Comparação dos resultados"
set xlabel "x"
set ylabel "φ(x)" rotate by 0
replot
set terminal wxt
set output


# with lines = w l
# linewidth = lw
# linecolor = lc
# dashtype = dt