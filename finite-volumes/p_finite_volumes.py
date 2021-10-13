import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
# Variáveis do domínio do problema e da simulação:
C = 0.8                       # Número de Courant
u = 0.5                       # Velocidade de propagação
N = 300                       # Número de nós na malha
L = 1.0                       # Domínio espacial
dx = L/N                      # Refinamento da discretização
dt = C*dx/u                   # Passo de tempo calculado
t0 = 0.0                      # Início da simulação
tf1 = 2.0                     # Tempo de interesse 1
tf2 = 5.0                     # Tempo de interesse 2
nsteps =(int)(tf1/dt)         # Número de passos de tempo
CFL = u*dt/dx

# Matrizes/Malhas:
X = np.linspace(0.0, L, N)  # Pontos em x para plotar
Q_up = np.zeros((nsteps, N))     # nsteps linhas (tempos) por N colunas (espaço)
Q_lax = np.zeros((nsteps, N))
Q_beam = np.zeros((nsteps, N))
Q_fromm = np.zeros((nsteps, N))

def show_parameters():
    print("\nNúmero de Courant: ", C)
    print("u: ", u)
    print("N: ", N)
    print("dx: ", dx)
    print("Tempo total: ", tf1)
    print("dt: ", dt)
    print("nsteps: ", nsteps)

def function_s(x):
    if (x>=0.6 and x<=0.8):
        return 1.0
    else:
        return 0.0

def function_phi(x):
    return math.exp(-200*((x - 0.3)**2)) + function_s(x)

def plot_x_por_y(x, y):
    plt.title("Condição Inicial")
    plt.plot(x, y, 'r', linewidth=2)
    plt.xlabel("x", fontsize = 11)
    plt.ylabel("Φ(x)", fontsize = 11)
    plt.grid(True, 'major', 'both')
    plt.show()

def plot_compare(x1, y1, x2, y2):
    plt.plot(x1, y1, 'r', label='Exato', linewidth=1)
    plt.plot(x2, y2, 'b', label='Numérico', linestyle='dashed',linewidth=2)
    plt.xlabel("x", fontsize = 11)
    plt.ylabel("Φ(x)", fontsize = 11)
    plt.legend()
    plt.grid(True, 'major', 'both')
    plt.show()
    
def plot_triplecompare(x, y1, y2, y3):
    plt.plot(x, y1, 'r', label='Inicial', linewidth=1)
    plt.plot(x, y2, 'b', label='Exato', linewidth=1)
    plt.plot(x, y3, linestyle='dashed', color='darkblue', label='Numérico')
    #plt.plot(x, y3, '', color='green', label='Numérico')
    plt.xlabel("x", fontsize = 11)
    plt.ylabel("Φ(x)", fontsize = 11)
    plt.legend()
    plt.grid(True, 'major', 'both')
    plt.savefig('Resultado.png')
    plt.show()
    
def solve_via_upwind(Q_up):
    # Itera no tempo:
    for n in range(1, nsteps):
        # Itera nas células espaciais:
        for i in range(0, N):
            # Periodicidade:
            if (i == 0):
                previous_pos = N-1
            else:
                previous_pos = i-1
            # Calcula os fluxos:
            Q_up[n, i] = Q_up[n-1, i] - C*(Q_up[n-1, i] - Q_up[n-1, previous_pos])

def solve_via_lax(Q_lax):
    # Itera no tempo:
    for n in range(1, nsteps):
        # Itera nas células espaciais:
        for i in range(0, N):
            # Periodicidade:
            if (i == 0):
                previous_pos = 0
            else:
                previous_pos = i-1
            if (i == N-1):
                next_pos = 0
            else:
                next_pos = i+1
            # Calcula os fluxos:
            Q_lax[n, i] = Q_lax[n-1, i] - 0.5*C*((Q_lax[n-1, next_pos] - Q_lax[n-1, previous_pos]) - C*(Q_lax[n-1, previous_pos] - 2*Q_lax[n-1, i] + Q_lax[n-1, next_pos]))
def solve_via_beam_warming(Q_beam):
    # Itera no tempo:
    for n in range(1, nsteps):
        # Itera nas células espaciais:
        for i in range(0, N):
            # Periodicidade:
            if (i == 0):
                previous_pos = 0
            else:
                previous_pos = i-1
            if (previous_pos == 0):
                p2_pos = 0
            else:
                p2_pos = previous_pos - 1
            # Calcula os fluxos:
            Q_beam[n, i] = Q_beam[n-1, i] - C*(Q_beam[n-1, i] - Q_beam[n-1, previous_pos]) - 0.5*C*(1-C)*(Q_beam[n-1, i] - 2*Q_beam[n-1, previous_pos] + Q_beam[n-1, p2_pos])
def solve_via_fromm(Q_fromm):
    # Itera no tempo:
    for n in range(1, nsteps):
        # Itera nas células espaciais:
        for i in range(0, N):
            # Periodicidade:
            if (i == N-1):
                next_pos = 0
            else:
                next_pos = i+1
            if (i == 0):
                previous_pos = 0
            else:
                previous_pos = i-1
            if (previous_pos == 0):
                p2_pos = 0
            else:
                p2_pos = previous_pos - 1
            # Calcula os fluxos:
            Q_fromm[n, i] = Q_fromm[n-1, i] - 0.25*C*((Q_fromm[n-1, next_pos] + 3*Q_fromm[n-1, i] - 5*Q_fromm[n-1, previous_pos] + Q_fromm[n-1, p2_pos]) - C*(Q_fromm[n-1, next_pos] - Q_fromm[n-1, i] - Q_fromm[n-1, previous_pos] + Q_fromm[n-1, p2_pos]))
# Preenche o tempo 0 com a condição inicial:
for x in range(0, N):
    Q_up[0, x] = function_phi(x*dx)
    Q_lax[0, x] = function_phi(x*dx)
    Q_beam[0, x] = function_phi(x*dx)
    Q_fromm[0, x] = function_phi(x*dx)

# Armazenaremos a condição incial, para comparação:
Ini = np.linspace(0.0, L, N)
for k in range(0, N):
    Ini[k] =  function_phi(k*dx)
# Solução exata no tempo tf:
Ext = np.linspace(0.0, L, N)
for k in range(0, N):
    Ext[k] =  function_phi(k*dx - u*(tf1-t0))

show_parameters()

solve_via_upwind(Q_up)
# solve_via_lax(Q_lax)
# solve_via_beam_warming(Q_beam)
# solve_via_fromm(Q_fromm)
plot_triplecompare(X, Ini, Ext, Q_up[nsteps-1])
# plot_triplecompare(X, Ini, Ext, Q_lax[nsteps-1])
# plot_triplecompare(X, Ini, Ext, Q_beam[nsteps-1])
# plot_triplecompare(X, Ini, Ext, Q_fromm[nsteps-1])


maxyvalue = (int)(np.amax(Q_up[nsteps-1]) )

timegif = plt.figure()
#plt.title("Evolucao da temperatura")
#plt.xlabel("Posição (m)", fontsize = 10)
#plt.ylabel("Temperatura (°C)", fontsize = 10)
ax = plt.axes(xlim=(0, L), ylim=(0, 1.1))
ax.grid()
time_text = ax.text(0.80, 0.03, '', transform=ax.transAxes)
line, = ax.plot([], [], 'b', lw=3)

def inicializacao():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text
def animacao(i):
    step = i*dt
    x = X
    y = Q_up[i]
    time_text.set_text('t = %.1f s' % step)
    line.set_data(x, y)
    return line, time_text

tsteps = (int)(nsteps/10)
anim = FuncAnimation(timegif, animacao, init_func=inicializacao, frames=nsteps, interval=10, blit=True)
anim.save('animacao.gif')
plt.show()

