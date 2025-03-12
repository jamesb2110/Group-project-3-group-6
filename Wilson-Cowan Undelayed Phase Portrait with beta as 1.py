import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from ddeint import ddeint  # Ensure you install ddeint via pip if not available

# Define the population activation function
def f(x, beta):  # Non-linear activation
    return 1 / (1 + np.exp(-beta * x))

# Parameters
pars = {
    'tau1': 0,    # Delay tau1
    'tau2': 0,    # Delay tau2
    'c1': -1,     # Connection strength for u -> u
    'c2': -0.4,      # Connection strength for u -> v
    'c3': -1,     # Connection strength for v -> u
    'c4': 0,      # Connection strength for v -> v
    'P': 0.65,         # External input for u
    'Q': 0.5,         # External input for v
    'beta': 1,    # Gain for the activation function
}

# Unpack parameters for convenience
c1, c2, c3, c4 = pars['c1'], pars['c2'], pars['c3'], pars['c4']
P, Q = pars['P'], pars['Q']
beta = pars['beta']
tau1, tau2 = pars['tau1'], pars['tau2']

# Define DDE system
def model(Y, t):
    # Y is the interpolated function (delayed states)
    u_t_tau1 = Y(t - tau1)[0]  # u(t - tau1)
    v_t_tau2 = Y(t - tau2)[1]  # v(t - tau2)
    u_t = Y(t)[0]              # u(t)
    v_t = Y(t)[1]              # v(t)
    
    du_dt = f(c1 * u_t_tau1 + c2 * v_t + P, beta) - u_t
    dv_dt = f(c3 * u_t + c4 * v_t_tau2 + Q, beta) - v_t
    
    return np.array([du_dt, dv_dt])

# Plot trajectories
plt.figure(figsize=(10, 10))

for a in np.linspace(1/7,6/7,5):
        # Define history function (initial conditions for t <= 0)
        def history_func(t):
            return np.array([a, 0])  # Initial state 

        # Time domain for simulation
        t = np.linspace(0, 20, 10000)  # Simulate from t=0 to t=20

        # Solve the delayed differential equations
        sol = ddeint(model, history_func, t)

        # Unpack solution
        u_t = sol[:, 0]  # u(t)
        v_t = sol[:, 1]  # v(t)
        plt.plot(u_t, v_t, color='black')
        
for a in np.linspace(1/7,6/7,5):
        # Define history function (initial conditions for t <= 0)
        def history_func(t):
            return np.array([a, 1])  # Initial state 

        # Time domain for simulation
        t = np.linspace(0, 20, 10000)  # Simulate from t=0 to t=20

        # Solve the delayed differential equations
        sol = ddeint(model, history_func, t)

        # Unpack solution
        u_t = sol[:, 0]  # u(t)
        v_t = sol[:, 1]  # v(t)
        plt.plot(u_t, v_t, color='black')
        
for a in np.linspace(1/7,6/7,5):
        # Define history function (initial conditions for t <= 0)
        def history_func(t):
            return np.array([0, a])  # Initial state 

        # Time domain for simulation
        t = np.linspace(0, 20, 10000)  # Simulate from t=0 to t=20

        # Solve the delayed differential equations
        sol = ddeint(model, history_func, t)

        # Unpack solution
        u_t = sol[:, 0]  # u(t)
        v_t = sol[:, 1]  # v(t)
        plt.plot(u_t, v_t, color='black')
        
for a in np.linspace(1/7,6/7,5):
        # Define history function (initial conditions for t <= 0)
        def history_func(t):
            return np.array([1, a])  # Initial state 

        # Time domain for simulation
        t = np.linspace(0, 20, 10000)  # Simulate from t=0 to t=20

        # Solve the delayed differential equations
        sol = ddeint(model, history_func, t)

        # Unpack solution
        u_t = sol[:, 0]  # u(t)
        v_t = sol[:, 1]  # v(t)
        plt.plot(u_t, v_t, color='black')

# Plot nullclines (reuse code from earlier for nullclines and vector field)
x = np.linspace(0, 1.1, 10)
y = np.linspace(0, 1.1, 10)
X, Y = np.meshgrid(x, y)  # Create the meshgrid for coordinates
u = f(c1 * (X - tau1) + c2 * Y + P, beta) - X  # X component with delay
v = f(c3 * X + c4 * (Y - tau2) + Q, beta) - Y  # Y component with delay

#plt.quiver(X, Y, u, v, color='black', alpha=0.5, pivot='middle',scale = 5)  # Vector field

# Plot nullclines
u_vals = np.linspace(0, 1.1, 1000)
v_vals = np.linspace(0, 1.1, 1000)
U, V = np.meshgrid(u_vals, v_vals)
U_nullcline = plt.contour(U, V, f(c1 * U + c2 * V + P, beta) - U, levels=[0], colors='b', linewidths=5, linestyles='dashed')
V_nullcline = plt.contour(U, V, f(c3 * U + c4 * V + Q, beta) - V, levels=[0], colors='magenta', linewidths=5, linestyles='dashed')

plt.plot(-1,-1,color='b', linewidth=3, linestyle='dashed', label='u-nullcline')
plt.plot(-1,-1,color='magenta', linewidth=3, linestyle='dashed', label='v-nullcline')
plt.plot(-1,-1,color='black', linewidth=3, label='Trajectories')
plt.scatter(0.4897970765124, 0.5025507087446,color='lime', label='Stable Node', zorder = 100, linewidth=3, edgecolors='black', s=400)

# Add labels and legend
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel('$u(t)$',fontsize=25)
plt.ylabel('$v(t)$',fontsize=25)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.savefig("Wilson-Cowan Undelayed Phase Portrait with beta as 1",dpi=500)
