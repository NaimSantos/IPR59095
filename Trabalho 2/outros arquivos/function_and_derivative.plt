set encoding utf8

#Variáveis
mu_o = 0.0008
mu_w = 0.001
C2 = mu_o/mu_w

# Funções e derivadas:
fw(x) = 1 / (1 + C2*(((1 - x)**2)/(x**2)))
gw(x)= (-(5*x-5)/(2*(x**3)))/ ((9*(x**2) - 10*x + 5)/(4*(x**2)))**2
pw(x) = (2*C2*((1-x)/x**2) + 2*C2*((1-x)**2)/(x**3)) / (1 + C2*(((1-x)**2)/(x**2)))**2


# Area de plotagem
set xrange[0:1]
plot fw(x) title "f_w(S_w)" w l lw 3 lc "blue", gw(x) title "f'_w(S_w)" w l lw 3 lc "red" #, pw(x) title "v2(S_w)" w l lw 3 lc "forest-green"
set terminal pngcairo enhanced
set output "Resultado.png"
set grid
#set title "Comparação dos resultados"
set xlabel "S_w"
set ylabel "F_w" rotate by 0
replot
set terminal wxt
set output
