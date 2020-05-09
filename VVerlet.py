import numpy as np
import matplotlib.pyplot as plt

#Setting functions for calculation of forces, potential energy,

def calculate_force(x, k):
    return -k*x

def calc_potential(x, k):
    return 0.5*k*x**2

def calc_kinetic(v, mass):
    return 0.5*mass*v**2

def vupdate_positions(x, v, dt, stepf=1.0):
    return x + v*dt*stepf

def vupdate_velocities(v, F, mass, dt, stepfraction=1.0):
    # type: (object, object, object, object, object) -> object
    return v + (0.5*dt*stepfraction)*F

#Function for Velocity Verlet MD

def VelocityVerletIntegrator(x, v, MDSteps=1000, dt=0.1, mass=1.0, k=10):

    x_trj = np.zeros(MDSteps, dtype=np.float64)
    v_trj = np.zeros(MDSteps, dtype=np.float64)

    x_trj[0] = x
    v_trj[0] = v

    for step in range(1, MDSteps):

        f = calculate_force(x, k)
        v = vupdate_velocities(v, f, mass, dt, stepfraction=0.5)
        x = vupdate_positions(x, v, dt, stepf=1.0)
        f = calculate_force(x, k)
        v = vupdate_velocities(v, f, mass, dt, stepfraction=0.5)
        x_trj[step] = x
        v_trj[step] = v

        x_trj = np.array(x_trj)
        v_trj = np.array(v_trj)
        Epot = calc_potential(x_trj, k)
        Ekin = calc_kinetic(v_trj, mass)
        Etot = Epot+Ekin

    return x_trj, v_trj, Epot, Ekin, Etot

x_0 = 1.0
v_0 = 1.0 #Saving initial velocities and coordinates to position and velocity vectors

MDSteps = int(1E05)
dt = 0.1

x, v, ep, ek, et = VelocityVerletIntegrator(x_0, v_0, MDSteps, dt, mass=1, k=10)

#Plot trajectory of phase space for 1D oscillator

fig = plt.figure(figsize=(12,6)) 

ax1 = fig.add_subplot(121)

ax1.plot(x, v)
ax1.set_xlabel("x")
ax1.set_ylabel("v")

ax2 = fig.add_subplot(122)

time = np.linspace(0, MDSteps, MDSteps*dt)

ep_short = [ep[i] for i in np.arange(0, len(ep), int(1/dt))]
ek_short = [ek[i] for i in np.arange(0, len(ek), int(1/dt))]
et_short = [et[i] for i in np.arange(0, len(et), int(1/dt))]

ax2.plot(time, ep_short, color="green", label="Potential Energy")
ax2.plot(time, ek_short, color="blue", label="Kinetic Energy")
ax2.plot(time, et_short, color="red", label="Total Energy")
ax2.legend(loc="best")

fig.tight_layout()

plt.show()
