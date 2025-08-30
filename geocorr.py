from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import math
from sympy import *

import numpy as np
from scipy.integrate import solve_ivp

#Definitions of functions and derivatives

def coth(a):
   return 1/(np.tanh(a))
def coshech(a):
   return 1/(np.sinh(a))

def u(y):
   return 1-2*y
def v(y):
   return 4*y*y

def f1(y,b):
    return (1/(b*b*y))*(coth(1/(v(y)))*coth(1/(v(y)))*(b*b*y*y -2) + coth(1/(v(y)))*(b*b*v(y)*v(y)-2*v(y))-8*b*b*y*y*y*y*y*(12*y+5)
             + b*b*v(y)*v(y)*(1+coth(1/(2*(y)))) + 2*y*y*y*b*b*(-3*coth(1/(2*y))+ coth(1/2*(y))*coth(1/(2*(y)))-1)
             -b*b*y*y*coth(1/(2*(y)))*coth(1/(2*(y)))
             + v(y)*v(y)*(v(y)-1) + 2*coth(1/(2*(y)))*(2*y+ coth(1/(2*(y))))+ y*y*b*b*coshech(1/v(y))*coshech(1/v(y)))


def my_ode_system(t, y):
    # Example: Null Geodesic
    # dy/dt = v
    # dv/dt = 3*y*y-y
    dydt = y[1]  # Velocity
    dvdt = 3*y[0]*y[0]-y[0] + 0.00000000001*f1(y[0],b) # Acceleration
    return [dydt, dvdt]

def my_ode_system1(t, y):
    dydt = y[1]  # Velocity
    dvdt = 3*y[0]*y[0]-y[0]  # Acceleration
    return [dydt, dvdt]


  
#interesting initial (1/4,1/0.4)(1/40,1/0.4)(1/3,1/30)
# Initial conditions: [initial position, initial velocity]
for i in range(50):
 y0 = [1/(3.0215),1/(370*(1.14214 + 0.00000015*(17.5)+0.000000015*i)) ]
 b=1/np.sqrt(y0[1]*y0[1]-2*y0[0]*y0[0]*y0[0]+y0[0]*y0[0])
 print(b)

# Time span for integration
 t_span = (0, 10*math.pi)

# Evaluate solution at specific time points for plotting
 t_eval = np.linspace(t_span[0], t_span[1], 100)

# Solve the ODE
 sol = solve_ivp(my_ode_system, t_span, y0, method='RK45', t_eval=t_eval)
 sol1 = solve_ivp(my_ode_system1, t_span, y0, method='RK45', t_eval=t_eval)
# Access the results
# solution.t contains the time points
# solution.y contains the corresponding y values (solution at each time point)
 r_positive = np.where(sol.y[0] > 0,1/sol.y[0], np.nan)
 r_corr=np.where(sol1.y[0] >0, 1/sol1.y[0], np.nan)

# Plot the solution
 plt.figure(figsize=(6,6)) 
 ax = plt.subplot(111, polar = True)
 ax.plot(sol.t,r_positive,c='green')
 ax.plot(sol1.t,r_corr,'--',c='red')

 plt.xlabel('Angle (phi)')
 plt.ylabel('Radius (r)')
 plt.grid(True)
 plt.show()
