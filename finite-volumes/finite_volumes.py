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
tf1 = 2.0                     # Tempo de interesse
nsteps =(int)(tf1/dt)         # Número de passos de tempo
CFL = u*dt/dx

# Matrizes/Malhas:
X = np.linspace(0.0, L, N)   # Pontos em x para plotar
Q_up = np.zeros((nsteps, N)) # nsteps linhas x N colunas
Q_lax = np.zeros((nsteps, N))
Q_beam = np.zeros((nsteps, N))
Q_fromm = np.zeros((nsteps, N))

Q_koren = np.zeros((nsteps, N))
Q_ospre = np.zeros((nsteps, N))
Q_van_albada = np.zeros((nsteps, N))

def function_s(x):
    if (x>=0.6 and x<=0.8):
        return 1.0
    else:
        return 0.0

def function_phi(x):
    return math.exp(-200*((x - 0.3)**2)) + function_s(x)

def phi_koren(teta):
    return max(0, min(2*teta, min((1 + 2*teta)/3, 2)))

def phi_ospre(teta):
    return 1.5*(teta**2 + teta)/(teta**2 + teta + 1)

def phi_albada(teta):
    return (teta**2 + teta)/(teta**2 + 1)

def function_phi(x):
    return math.exp(-200*((x - 0.3)**2)) + function_s(x)
    
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

def solve_with_minmod(Q):
    # Itera no tempo:
    for n in range(1, nsteps):
        # Itera nas células espaciais:
        for i in range(0, N):
            # Periodicidade:
            if (i == 0):
                i_prev = N-1
            else:
                i_prev = i-1
            if (i_prev == 0):
                i_prev_2 = N-1
            else:
                i_prev_2 = i_prev-1
            if (i == N-1):
                i_next = 0
            else:
                i_next = i+1
            if (i_next == N-1):
                i_next_2 = 0
            else:
                i_next_2 = i_next+1

            if (u > 0):
                teta_prev = (Q[n-1, i_prev] - Q[n-1, i_prev_2])/(Q[n-1, i] - Q[n-1, i_prev])
                teta_next = (Q[n-1, i] - Q[n-1, i_prev])/(Q[n-1, i_next] - Q[n-1, i])
                Q[n, i] = Q[n-1, i] - C*(Q[n-1, i] - Q[n-1, i_prev]) - 0.5*C*(1 - C)*(phi_koren(teta_next)*(Q[n-1, i_next] - Q[n-1, i]) - phi_koren(teta_prev)*(Q[n-1, i] - Q[n-1, i_prev]))
            else:
                teta_prev = (Q[n-1, i_next] - Q[n-1, i])/(Q[n-1, i] - Q[n-1, i_prev])
                teta_next = (Q[n-1, i_next_2] - Q[n-1, i_next])/(Q[n-1, i_next] - Q[n-1, i])
                Q[n, i] = Q[n-1, i] - C*(Q[n-1, i_next] - Q[n-1, i]) + 0.5*C*(1 + C)*(phi_koren(teta_next)*(Q[n-1, i_next] - Q[n-1, i]) - phi_koren(teta_prev)*(Q[n-1, i] - Q[n-1, i_prev]))

def solve_with_superbee(Q):
    # Itera no tempo:
    for n in range(1, nsteps):
        # Itera nas células espaciais:
        for i in range(0, N):
            # Periodicidade:
            if (i == 0):
                i_prev = N-1
            else:
                i_prev = i-1
            if (i_prev == 0):
                i_prev_2 = N-1
            else:
                i_prev_2 = i_prev-1
            if (i == N-1):
                i_next = 0
            else:
                i_next = i+1
            if (i_next == N-1):
                i_next_2 = 0
            else:
                i_next_2 = i_next+1

            if (u > 0):
                teta_prev = (Q[n-1, i_prev] - Q[n-1, i_prev_2])/(Q[n-1, i] - Q[n-1, i_prev])
                teta_next = (Q[n-1, i] - Q[n-1, i_prev])/(Q[n-1, i_next] - Q[n-1, i])
                Q[n, i] = Q[n-1, i] - C*(Q[n-1, i] - Q[n-1, i_prev]) - 0.5*C*(1 - C)*(phi_ospre(teta_next)*(Q[n-1, i_next] - Q[n-1, i]) - phi_ospre(teta_prev)*(Q[n-1, i] - Q[n-1, i_prev]))
            else:
                teta_prev = (Q[n-1, i_next] - Q[n-1, i])/(Q[n-1, i] - Q[n-1, i_prev])
                teta_next = (Q[n-1, i_next_2] - Q[n-1, i_next])/(Q[n-1, i_next] - Q[n-1, i])
                Q[n, i] = Q[n-1, i] - C*(Q[n-1, i_next] - Q[n-1, i]) + 0.5*C*(1 + C)*(phi_ospre(teta_next)*(Q[n-1, i_next] - Q[n-1, i]) - phi_ospre(teta_prev)*(Q[n-1, i] - Q[n-1, i_prev]))
            
def show_parameters():
    print("\nNúmero de Courant: ", C)
    print("u: ", u)
    print("N: ", N)
    print("dx: ", dx)
    print("Tempo total: ", tf1)
    print("dt: ", dt)
    print("nsteps: ", nsteps)
            
def plot_x_por_y(x, y):
    plt.title("Condição Inicial")
    plt.plot(x, y, 'r', linewidth=2)
    plt.xlabel("x", fontsize = 11)
    plt.ylabel("Φ(x)", fontsize = 11)
    plt.grid(True, 'major', 'both')
    plt.show()

def plot_compare(x, y1, y2, file_name):
    plt.plot(x, y1, 'r', label='Solução Exata', linewidth=2)
    plt.plot(x, y2, 'b', label=file_name, linestyle='dashed',linewidth=2)
    plt.xlabel("x", fontsize = 11)
    plt.ylabel("Φ(x)", fontsize = 11)
    plt.legend()
    plt.savefig(file_name + '.png')
    plt.grid(True, 'major', 'both')
    plt.show()
    
def plot_triplecompare(x, y1, y2, y3):
    plt.plot(x, y1, 'r', label='Inicial', linewidth=1)
    plt.plot(x, y2, 'b', label='Exato', linewidth=1)
    plt.plot(x, y3, linestyle='dashed', color='darkblue', label='Numérico')
    plt.xlabel("x", fontsize = 11)
    plt.ylabel("Φ(x)", fontsize = 11)
    plt.legend()
    plt.grid(True, 'major', 'both')
    plt.savefig('Resultado.png')
    plt.show()

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

# solve_with_superbee(Q_ospre)
# plot_compare(X, Ini, Q_ospre[nsteps-1], 'Método Ospre')
solve_via_upwind(Q_up)
plot_triplecompare(X, Ini, Ext, Q_up[nsteps-1])


################# Animação #################
# timegif = plt.figure()
# ax = plt.axes(xlim=(0, L), ylim=(-0.3, 1.3))
# ax.grid()
# time_text = ax.text(0.85, 0.95, '', transform=ax.transAxes)
# line, = ax.plot([], [], 'b', lw=3)

# def inicializacao():
    # line.set_data([], [])
    # time_text.set_text('')
    # return line, time_text
# def animacao(i):
    # step = i*dt
    # x = X
    # y = Q_up[i]
    # time_text.set_text('t = %.1f s' % step)
    # line.set_data(x, y)
    # ax.set_title("Evolução da solução  (Upwind)")
    # ax.axes.set_xlabel("x", fontsize = 10)
    # ax.axes.set_ylabel("Φ(x)", fontsize = 10)
    # return line, time_text

# tsteps = (int)(nsteps/10)
# anim = FuncAnimation(timegif, animacao, init_func=inicializacao, frames=nsteps, interval=20, blit=True)
# anim.save('Animacao1.mp4')
# plt.show()
