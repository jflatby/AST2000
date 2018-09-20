import numpy as np
from ast2000solarsystem_27_v4 import AST2000SolarSystem
import matplotlib.pyplot as plt
from scipy.interpolate import spline



solar_system = AST2000SolarSystem(11466)

planet_mass = solar_system.mass[6]                                  #[Solar masses]
planet_radius = solar_system.radius[6]
semi_major_axes = solar_system.a[6]

gravitational_constant = 4 * np.pi**2                               #[AU^3 * Yr^(-2) * solar masses^(-1)]
star_mass = solar_system.star_mass
star_radius = solar_system.star_radius
reduced_mass = star_mass*planet_mass/(star_mass + planet_mass)

planet_orbital_period = np.sqrt((semi_major_axes**3) / (star_mass + planet_mass))
total_time = planet_orbital_period*5
print total_time
time_steps_per_year = 10000
total_time_steps = int(time_steps_per_year * total_time)
dt = total_time/total_time_steps

positions = np.zeros((2, 2, total_time_steps))
velocities = np.zeros((2, 2, total_time_steps))
accelerations = np.zeros((2, 2, total_time_steps))

positions[:, 0, 0] = [0, 0]
positions[:, 1, 0] = [solar_system.x0[6], solar_system.y0[6]]

velocities[:, 1, 0] = [solar_system.vx0[6], solar_system.vy0[6]]

mass_center = np.zeros((2, total_time_steps))
system_energy = np.zeros(total_time_steps)

accelerations[:, 0, 0] = -gravitational_constant*planet_mass *(positions[:, 0, 0] - positions[:, 1, 0]) / np.linalg.norm((positions[:, 0, 0] - positions[:, 1, 0]))**3
accelerations[:, 1, 0] = -gravitational_constant*star_mass *(positions[:, 1, 0] - positions[:, 0, 0]) / np.linalg.norm((positions[:, 1, 0] - positions[:, 0, 0]))**3

mass_center[:, 0] = star_mass/(star_mass + planet_mass) * positions[:, 0, 0] + planet_mass/(star_mass + planet_mass) * positions[:, 1, 0]
system_energy[0] = 0.5 * reduced_mass * \
                       np.linalg.norm(velocities[:, 0, 0] + velocities[:, 1, 0]) ** 2 \
                       - gravitational_constant * (star_mass + planet_mass) * \
                         reduced_mass / np.linalg.norm(mass_center[:, 0])

for t in range(total_time_steps-1):

    r_vector = positions[:, 0, t] - positions[:, 1, t]
    r = np.linalg.norm(r_vector)

    #Star
    positions[:, 0, t+1] = positions[:, 0, t] + velocities[:, 0, t]*dt + 0.5*accelerations[:, 0, t]*dt**2
    accelerations[:, 0, t+1] = -gravitational_constant * planet_mass * r_vector / r**3
    velocities[:, 0, t+1] = velocities[:, 0, t] + 0.5*(accelerations[:, 0, t] + accelerations[:, 0, t+1]) * dt

    #Planet
    positions[:, 1, t+1] = positions[:, 1, t] + velocities[:, 1, t] * dt + 0.5 * accelerations[:, 1, t] * dt ** 2
    accelerations[:, 1, t+1] = -gravitational_constant * star_mass * -r_vector / r**3
    velocities[:, 1, t+1] = velocities[:, 1, t] + 0.5*(accelerations[:, 1, t] + accelerations[:, 1, t + 1]) * dt

    #Energy sdjfnlk
    mass_center[:, t+1] = star_mass/(star_mass + planet_mass) * positions[:, 0, t+1] + planet_mass/(star_mass + planet_mass) * positions[:, 1, t+1]

    system_energy[t+1] = 0.5 * reduced_mass * \
                       np.linalg.norm(velocities[:, 0, t] + velocities[:, 1, t]) ** 2 \
                       - gravitational_constant * (star_mass + planet_mass) * \
                         reduced_mass / np.linalg.norm(mass_center[:, t])


print system_energy

#plt.plot(positions[0, 0, :], positions[1, 0, :], positions[0, 1, :], positions[1, 1, :], mass_center[0, :], mass_center[1, :])
#plt.show()

t = np.linspace(0, total_time, total_time_steps)

plt.style.use("ggplot")

plt.xlabel("Earth years")
plt.ylabel("AU / Yr")

i = np.pi/4
v_peculiar = 2        #AU/yr
obs = np.pi/4

uvec = "teit" # Sikkert feil np.array([np.cos(i) * np.sin(obs), np.sin(i), np.cos(i) * np.cos(obs)])
v_rstar = np.einsum('ij,j->i',velocities[:, 0, :], uvec)

plt.plot(t, v_rstar)
plt.show()

v_r_max = np.amax(velocities[0, 0, :])

y_noise = velocities[0, 0, :] + np.random.normal(0, 0.2*v_r_max, total_time_steps)

plt.plot(t, y_noise, linewidth=0.1)
plt.show()


# sjekker nar planeten er foran sola med x = uendelig og y = 0 som ref. pkt
"""

x = np.linspace(2000, 2300, 300)
y = [1 for i in range(300)]

for t in range(2000, 2300):
    if positions[0, 1, t] > positions[0, 0, t]: #If planet is on this side of star
        if positions[1, 1, t] - planet_radius*6.68459e-9 <= positions[1, 0, t] + star_radius*6.68459e-9 and positions[1, 1, t] + planet_radius*6.68459e-9 >= positions[1, 0, t] - star_radius*6.68459e-9:
            if positions[1, 1, t] + planet_radius*6.68459e-9 <= positions[1, 0, t] + star_radius*6.68459e-9 and positions[1, 1, t] - planet_radius*6.68459e-9 >= positions[1, 0, t] - star_radius*6.68459e-9:
                y[t-2000] = 1 - (planet_radius*6.68459e-9)**2 / (star_radius*6.68459e-9)**2
                #print "Innenfor"
            else:
                print "Delvis innenfor"
"""

