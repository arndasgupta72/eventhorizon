import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Defining the physical constants

G = 6.67*10**-11                    # Universal Gravitational constant
M = 2*10**30                        # Mass of Sun
AU = 1.496*10**11                   # Astronomical Unit (A.U.)
c = 299792458                       # Speed of light
r_s = 2*G*M/c**2                    # Schwarzschild radius
L_earth = 4.4536*10**15             # Angular momentum of Earth per unit mass
L_mercury = 2.7701*10**15           # Angular momentum of Mercury per unit mass
merc_per = 0.307499*AU              # Perihelion distance of Mercury
earth_per = 0.9832899*AU            # Perihelion distance of Earth
e_merc = 0.21                       # eccentricity of mercury's orbit (in Newtonian mechanics)
e_earth = 0.017                     # eccentricity of earth's orbit (in Newtonian mechanics)

# Defining the derivative function
# here x = angular position (phi)
# and y = radial position (r)

# Defining the derivative function
# here x = angular position (phi)
# and y = radial position (r)

def ddx(x,y):
    u = y[0]
    v = y[1]
    return np.array([v,3*G*M*u**2/c**2-u])

# Defining RK4 method

def rk4(f,x,y,h):
    k1 = h*f(x,y)
    k2 = h*f(x+h/2,y+k1/2)
    k3 = h*f(x+h/2,y+k2/2)
    k4 = h*f(x+h,y+k3)
    return y+(k1+2*k2+2*k3+k4)/6

# Defining the Caller function

def Caller(mtd,f,xs,y_ini,h):
    N = len(xs)
    y = y_ini
    ys = np.zeros((N,2),float)
    for i in range(N):
        x = xs[i]
        ys[i] = y
        y = mtd(ddx,x,y,h)
    return ys

# Initial Conditions

# Initial Coditions

r0 = 1.50001*r_s                              # Initial Radial Position
u0 = 1/r0
v0 = 0                                          # Initial Radial Velocity
y_ini = np.array([u0,v0])

# Time step

h = 0.0001

# Getting a solution

n = 2.0                                        # Number of full circles

xs = np.arange(0,2*np.pi*n,h)                # Stores the angular coordinates
ys = Caller(rk4,ddx,xs,y_ini,h)              # Stores the radial coordinates and radial velocities

r = 1/ys[:,0]                                  # Stores just the radial coordinates
r = r/(1*AU)                                 # Normalizes the radial distance w.r.t. 1 A.U.

# Compute x and y coordinates
x_ax = r * np.cos(xs)
y_ax = r * np.sin(xs)

# Animation
fig, ax = plt.subplots()

# Create a point (ball) to represent the planet
planet, = ax.plot([], [], 'bo', ms=5)  # 'bo' means blue color and 'ms' stands for marker size

# Plot Event Horizon and additional elements

pps = np.arange(0,2*np.pi+0.01,0.01)
xxs = [r_s/AU*np.cos(i) for i in pps]
yys = [r_s/AU*np.sin(i) for i in pps]

l = 10**-7.2

plt.plot(xxs,yys,c='red',linestyle = 'dotted',label="Schwarzschild Radius")
plt.scatter([0], [0], marker='.', c='black', label="Sun")
plt.plot(x_ax,y_ax,c='orange',linestyle = 'dashed',label="Photon's Path")
plt.xlabel("Distance (in A.U.)")
plt.ylabel("Distance (in A.U.)")
plt.title("Photon Orbit in Schwarzschild Metric")
plt.grid()
plt.xlim(-l,l)
plt.ylim(-l,l)
plt.legend()

def init():
    planet.set_data([], [])
    return planet,

def update(frame):
    # Update the position of the planet
    planet.set_data([x_ax[frame]], [y_ax[frame]])
    return planet,

animation = FuncAnimation(fig, update, frames=np.arange(0, len(x_ax), 700), init_func=init, interval=10, blit=True)

# Save the animation

# writer = PillowWriter(fps=20)
# animation.save("Animation_1.gif", writer=writer)

plt.show()