set key right bottom

set style line 1 linecolor "red" linetype 1 linewidth 2 pointtype 1 pointsize 2
set style line 2 linecolor "magenta" linetype 2 linewidth 2 pointtype 2 pointsize 2
set style line 3 linecolor "blue" linetype 1 linewidth 2 pointtype 8 pointsize 2
set style line 4 linecolor "dark-green" linetype 1 linewidth 2 pointtype 6 pointsize 2
set style line 5 linecolor "gold" linetype 1 linewidth 2 pointtype 3 pointsize 2


plot "dat_phi_015.txt" using 1:2 with lines ls 1 title "Φ = 0,15", "dat_phi_020.txt" using 1:2 with lines ls 2 title "Φ = 0,20", "data_padrao.txt" using 1:2 with lines ls 3 title "Caso Padrão", "dat_phi_030.txt" using 1:2 with lines ls 4 title "Φ = 0,30", "dat_phi_035.txt" using 1:2 with lines ls 5 title "Φ = 0,35"

set terminal pngcairo enhanced
set output "pressao_no_espaco.png"
set grid xtics ytics lw 2
set grid
set xlabel "x (m)"
set ylabel "Pressao (kPa)"
replot
set terminal wxt
set output

# lc - linecolor
# lt - linetype
# lw - linewidth
# pt - pointtype
# ps - pointsize
# w  - with
# lp - linespoints
# ls - linestyle
# w l -> with lines
# lw ->  linewidth
# lc -> linecolor 
# dt -> dashtype
# rotate by 0 -> para cada eixo (xlabel, ylabel), rotaciona o título do eixo
# set key outside: legenda fora da área plotável. 



# plot "data_padrao.txt" using 1:2 w l lw 2 lc "blue" title "Caso Padrão", "data_dt1.txt" using 1:2 w l lw 2 lc "gold" title " dt = 1,0 dia(s)",  "data_dt2.5.txt" using 1:2 w l lw 2 lc "red" title "dt = 2,5 dia(s)", "data_dt5.txt" using 1:2 w l lw 2 lc "dark-green" title "dt = 5,0 dia(s)"