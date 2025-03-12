import numpy as np
import matplotlib.pyplot as plt
from ddeint import ddeint

# Define parameters
N = 2
tau_1 = 0.7
tau_2 = 0.9
rho = 2.4
epsilon = -0.3
P = 0.65
Q = 0.5
c1 = -1.0
c2 = -0.4
c3 = -1.0
c4 = 0.0
beta = 60
u_initial = [0.1,0.9] #np.linspace(1/(N+1),1-1/(N+1),N)
v_initial = [0.9,0.1] #np.linspace(1/(N+1),1-1/(N+1),N)
color_array = ["blue","magenta"]

# Define the logistic function f(x)
def f(x, beta):
    return 1 / (1 + np.exp(-beta * x))

# Define the system of DDEs as a function
def wilson_cowan_ring(Y, t, rho, epsilon, P, Q, c1, c2, c3, c4, W, beta, N, tau_1, tau_2):
    u = Y(t)[:N]  # excitatory population values (u_i)
    v = Y(t)[N:]  # inhibitory population values (v_i)

    # Delayed values
    u_delayed = Y(t - tau_1)[:N]
    v_delayed = Y(t - tau_2)[N:]

    # Initialize derivatives
    du_dt = np.zeros(N)
    dv_dt = np.zeros(N)

    for i in range(N):
        # Sum over all j for the coupling terms
        coupling_sum = sum(W[i, j] * u_delayed[j] for j in range(N))

        du_dt[i] = -u[i] + f(c1 * u_delayed[i] + c2 * v_delayed[i] + P + epsilon * coupling_sum, beta)
        dv_dt[i] = -v[i] + f(c3 * u_delayed[i] + c4 * v_delayed[i] + Q, beta)

    return np.concatenate([du_dt, dv_dt])

# Define connectivity matrices
def get_matrix(N, matrix_type):
    W = np.zeros((N, N))
    if matrix_type == "RING":
        for i in range(0,N-1):
            W[i,i+1] = 0.5
            W[i,i-1] = 0.5
        W[N-1,0] = 0.5
        W[N-1,N-2] = 0.5
    elif matrix_type == "STAR":
        W.fill(1/(N-1))
        np.fill_diagonal(W, 0)
    elif matrix_type == "CENTRE":
        W[0,1:] = 1/(N-1)
        W[1:,0] = 1
    return W

# Function to solve the system using ddeint
def solve_system(u_initial, v_initial, N, tau_1, tau_2, rho, epsilon, P, Q, c1, c2, c3, c4, W, beta, tmin=0.0, tmax=50.0, timesteps=1000):
    t = np.linspace(tmin, tmax, timesteps)

    def history(t):       
        return np.concatenate([u_initial, v_initial])

    system = lambda y, t: wilson_cowan_ring(y, t, rho, epsilon, P, Q, c1, c2, c3, c4, W, beta, N, tau_1, tau_2)
    
    # Solve using ddeint
    sol = ddeint(system, history, t)
    return t, sol

# Main simulation function
def run_simulation(N, tau_1, tau_2, rho, epsilon, P, Q, c1, c2, c3, c4, beta, matrix_type, u_initial, v_initial):
    W = get_matrix(N, matrix_type)
    
    tmin, tmax, timesteps = 0.0, 50.0, 1000

    t, sol = solve_system(u_initial, v_initial, N, tau_1, tau_2, rho, epsilon, P, Q, c1, c2, c3, c4, W, beta, tmin, tmax, timesteps)
        
    u_sol = sol[:, :N]
    v_sol = sol[:, N:]

    plt.figure(num=1,figsize=(10, 10))
    for i in range(N):
        plt.plot(u_sol[:, i], v_sol[:, i], label=f"$u_{i+1}, v_{i+1}$",color=color_array[i],linewidth=3)
        plt.scatter(u_initial[i],v_initial[i],color=color_array[i], zorder = 100, linewidth=3, edgecolors='black', s=200)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('$u(t)$',fontsize=25)
    plt.ylabel('$v(t)$',fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.savefig("Wilson-Cowan 2-Node Network Phase Portrait",dpi=500)
    
    plt.figure(num=2,figsize=(20, 10))
    for i in range(N):
        plt.plot(t, u_sol[:, i], label=f"$u_{i+1}(t)$", linewidth=4,color="blue")
        plt.plot(t, v_sol[:, i], label=f"$v_{i+1}(t)$", linewidth=4,color="magenta")
    plt.xlim(0, 20)
    plt.ylim(0, 1)
    plt.xlabel('$t$',fontsize=25)
    plt.ylabel('$u(t)/v(t)$',fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.savefig("Wilson-Cowan 2-Node Network u_i(t) and v_i(t)",dpi=500)

    plt.show()

# Example usage

run_simulation(N, tau_1, tau_2, rho, epsilon, P, Q, c1, c2, c3, c4, beta, "CENTRE", u_initial, v_initial)
