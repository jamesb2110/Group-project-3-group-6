import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
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

pre_bifurcation_u = np.array([])
pre_bifurcation_v = np.array([])
post_bifurcation_u = np.array([])
post_bifurcation_v = np.array([])
unstable_u = np.array([])
unstable_v = np.array([])

# Unpack parameters for convenience
c1, c2, c3, c4 = pars['c1'], pars['c2'], pars['c3'], pars['c4']
P, Q = pars['P'], pars['Q']
beta = pars['beta']
tau1, tau2 = pars['tau1'], pars['tau2']

for beta in np.linspace(22.69375,60,100):
    def system(uv):
        return np.array([f(c1 * uv[0] + c2 * uv[1] + P, beta) - uv[0],f(c3 * uv[0] + c4 * uv[1] + Q, beta) - uv[1]])
    unstable_u = np.append(unstable_u,fsolve(system,np.array([0.5,0.5]))[0])
    unstable_v = np.append(unstable_v,fsolve(system,np.array([0.5,0.5]))[1])

for beta in np.linspace(1,60,250):
    
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
    
    # Define history function (initial conditions for t <= 0)
    def history_func(t):
        return np.array([0.5, 1])  # Initial state 
    
    # Time domain for simulation
    t = np.linspace(0, 20, 1000)  # Simulate from t=0 to t=20
   
    # Solve the delayed differential equations
    sol = ddeint(model, history_func, t)
    
    # Unpack solution
    pre_bifurcation_u = np.append(pre_bifurcation_u,sol[-1, 0])
    pre_bifurcation_v = np.append(pre_bifurcation_v,sol[-1,1])
    
for beta in np.linspace(22.69375,60,100):
    
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
    
    # Define history function (initial conditions for t <= 0)
    def history_func(t):
        return np.array([0.5, 0])  # Initial state 
    
    # Time domain for simulation
    t = np.linspace(0, 20, 1000)  # Simulate from t=0 to t=20
   
    # Solve the delayed differential equations
    sol = ddeint(model, history_func, t)
    
    # Unpack solution
    post_bifurcation_u = np.append(post_bifurcation_u,sol[-1, 0])
    post_bifurcation_v = np.append(post_bifurcation_v,sol[-1,1])

plt.figure(figsize=(10, 10))
plt.plot(np.linspace(1,60,250),pre_bifurcation_u,color='black',linewidth=5)
plt.plot(np.linspace(22.69375,60,100),post_bifurcation_u,color='black',linewidth=5)
plt.plot(np.linspace(22.69375,60,100),unstable_u,color='black',linewidth=5,linestyle='dashed')
plt.plot(np.array([22.69375,22.69375]),np.array([0,1]),color='r',linestyle='dashed',linewidth=5)
plt.xlim(0, 60)
plt.ylim(0, 1)
plt.xlabel(r'$\beta$',fontsize=25)
plt.ylabel('$u$',fontsize=25)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.savefig("Un-delayed bifurcation diagram (u values)",dpi=500)