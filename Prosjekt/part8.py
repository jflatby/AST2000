from ast2000solarsystem_27_v5 import AST2000SolarSystem
import numpy as np
import matplotlib.pyplot as plt

star_system = AST2000SolarSystem(11466)

#star_system.part2B_4(1)

c = 1
L0 = 200
v0 = 0.99
g = 3.34e-10
gamma = np.sqrt(1-v0**2)
deltaT = L0/v0 - v0/g

total_time = deltaT
number_of_time_steps = 1000
dt = total_time / number_of_time_steps

ty = np.linspace(0, L0/v0 + deltaT, number_of_time_steps)
ty_ = np.zeros(number_of_time_steps)

v = v0

for t in xrange(number_of_time_steps):
    if t*dt < L0/v0:
        ty_[t] = ty[t]/np.sqrt(1-v0**2/c**2)
        print ty[t], ty_[t]
    else:
        v = v + g*dt
        ty_[t] = -L0 * v + L0 / v + (t*dt - L0/v )

#plt.axes().set_aspect('equal', 'datalim')

plt.plot(ty, ty_)
plt.show()