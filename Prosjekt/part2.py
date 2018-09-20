from __future__ import division
import matplotlib.pyplot as plt
from ast2000solarsystem_27_v4 import AST2000SolarSystem
import numpy as np
plt.style.use("ggplot")

solar_system = AST2000SolarSystem(11466)

gravitational_constant = 4 * np.pi**2                               #[AU^3 * Yr^(-2) * solar masses^(-1)]
number_of_planets = solar_system.number_of_planets
semi_major_axes = solar_system.a                                    #[AU]
star_mass = solar_system.star_mass                                  #[solar masses]
star_radius = solar_system.star_radius                              #[km]
planet_masses = solar_system.mass                                   #[solar masses]
planet_radii = solar_system.radius                                  #[km]

print solar_system.e

home_planet_orbital_period = np.sqrt((semi_major_axes[0]**3) / (star_mass + planet_masses[0]))
#Planet orbital period: np.sqrt((semi_major_axes[p]**3) / (star_mass + planet_masses[0]))
total_time = np.sqrt((semi_major_axes[5]**3) / (star_mass + planet_masses[5])) #Years
print total_time
time_steps_per_year = 70000
total_time_steps = int(time_steps_per_year*total_time)
dt = total_time/total_time_steps

times = np.linspace(0, total_time, total_time_steps)

planet_positions = np.zeros((2, number_of_planets, total_time_steps))         #[AU]
planet_velocities = np.zeros((2, number_of_planets, total_time_steps))        #[AU/Yr]
planet_accelerations = np.zeros((2, number_of_planets, total_time_steps))
distance_to_star = np.zeros(number_of_planets)

planet_apoapsis = np.zeros((number_of_planets, 2))
planet_periapsis = np.zeros((number_of_planets, 2)) + 10000

au_to_m = 1.496e+11

for p in range(number_of_planets):
    planet_positions[:, p, 0] = [solar_system.x0[p], solar_system.y0[p]]     #[AU]
    planet_velocities[:, p, 0] = [solar_system.vx0[p], solar_system.vy0[p]]  #[AU/Yr]
    distance_to_star[p] = np.sqrt(planet_positions[0, p, 0]**2 + planet_positions[1, p, 0]**2)
    planet_accelerations[:, p, 0] = -gravitational_constant*star_mass * planet_positions[:, 0, 0]/distance_to_star[p]**3


for t in range(total_time_steps-1):
        planet_positions[:, :, t+1] = planet_positions[:, :, t] + planet_velocities[:, :, t]*dt + 0.5*planet_accelerations[:, :, t]*dt**2
        r = np.sqrt(planet_positions[0, :, t] ** 2 + planet_positions[1, :, t] ** 2)
        planet_accelerations[:, :, t+1] = -gravitational_constant * star_mass * planet_positions[:, :, t+1] / r ** 3
        planet_velocities[:, :, t+1] = planet_velocities[:, :, t] + 0.5*(planet_accelerations[:, :, t] + planet_accelerations[:, :, t+1])*dt

        #Check if at heighest or lowest point, and update array accordingly:
        for p in range(number_of_planets):
            if np.linalg.norm(planet_positions[:, p, t]) > planet_apoapsis[p, 1]:
                planet_apoapsis[p] = [t, np.linalg.norm(planet_positions[:, p, t])]
            if np.linalg.norm(planet_positions[:, p, t]) < planet_periapsis[p, 1]:
                planet_periapsis[p] = [t, np.linalg.norm(planet_positions[:, p, t])]



        if t%70000 == 0:
            print t/70000

planet_apoapsis = np.array(
    [[  0.00000000e+00,   4.44345620e+00],
     [  1.80242000e+05,   6.39799144e+00],
     [  7.27839000e+05,   2.74951134e+01],
     [  3.69806200e+06,   3.96877237e+01],
     [  2.55220000e+04,   1.94768310e+01],
     [  5.27633000e+05,   1.50277315e+01],
     [  5.64450000e+04,   2.84413133e+00]])

planet_periapsis = np.array(
    [[  2.03413000e+05,   4.26463310e+00],
     [  5.33625000e+05,   6.18645634e+00],
     [  3.95221800e+06,   2.74548341e+01],
     [  9.07585500e+06,   3.75932729e+01],
     [  1.89583700e+06,   1.87422095e+01],
     [  1.71915000e+06,   1.32690266e+01],
     [  1.61402000e+05,   2.75790904e+00]]
)

def analytical_orbit(planet):
    analytical_position = np.zeros(total_time_steps)
    analytical_position[0] = np.sqrt(solar_system.x0[planet]**2 + solar_system.y0[planet]**2)
    theta = np.linspace(0, 2*np.pi, total_time_steps)

    a = solar_system.a[planet]
    e = solar_system.e[planet]
    f = theta - solar_system.psi[planet]

    for t in range(total_time_steps-1):
        analytical_position[t+1] = a*(1-e**2)/(1+e*np.cos(f[t+1]))

    return analytical_position*np.cos(theta), analytical_position*np.sin(theta)

def get_angle_at_time(planet, t):
    pos_vec = planet_positions[:, planet, t]
    print pos_vec
    return np.arctan2(pos_vec[1], pos_vec[0])


def integrate_area(p):
    N = 100
    t_apoapsis = planet_apoapsis[p, 0]
    T = 1
    t_0 = int(t_apoapsis/70000 - T/2)
    t_n = int(t_apoapsis/70000 + T/2)
    dt = T / N
    theta = np.linspace(get_angle_at_time(p, t_0), get_angle_at_time(p, t_n), N)
    dAdt = 0


    dtheta = theta[0] - theta[-1]
    w = dtheta / T

    for t in range(N):
        r = np.linalg.norm(planet_positions[:, p, int(t*dt)])
        v_t = theta[t]/dt
        dAdt += 0.5*r**2 * v_t
    print "Apoapsis: dA/dt = ", dAdt

    t_periapsis = planet_periapsis[p, 0]
    t_0 = int(t_periapsis/70000 - T/2)
    t_n = int(t_periapsis/70000 + T/2)
    theta = np.linspace(get_angle_at_time(p, t_0), get_angle_at_time(p, t_n), N)

    dAdt = 0
    for t in range(N):
        r = np.linalg.norm(planet_positions[:, p, int(t*dt)])
        v_t = theta[t]/dt
        dAdt += 0.5*r**2 * v_t
    print "Periapsis: dA/dt = ", dAdt

    ##Kepler integration:
    t_apoapsis = planet_apoapsis[p, 0]
    t_periapsis = planet_periapsis[p, 0]
    t_0_a = int(t_apoapsis - 0.2 * 70000)
    t_n_a = int(t_apoapsis + 0.2 * 70000)
    t_0_p = int(t_periapsis - 0.2 * 70000)
    t_n_p = int(t_periapsis + 0.2 * 70000)

    plt.figure()
    apo = plt.Polygon(np.array([planet_positions[:, p, t_0_a], planet_positions[:, p, t_n_a], np.array([0, 0])]),
                      color="blue")
    plt.gca().add_patch(apo)
    peri = plt.Polygon([planet_positions[:, p, t_0_p], planet_positions[:, p, t_n_p], [0, 0]], color="red")
    plt.gca().add_patch(peri)


#integrate_area(0)

#solar_system.check_planet_positions(planet_positions, total_time, time_steps_per_year)
#solar_system.orbit_xml(planet_positions, times)


plt.plot(planet_positions[0, 0, :], planet_positions[1, 0, :])
#plt.plot(planet_positions[0, 1, :], planet_positions[1, 1, :])
#plt.plot(planet_positions[0, 2, :], planet_positions[1, 2, :])
#plt.plot(planet_positions[0, 3, :], planet_positions[1, 3, :])
#plt.plot(planet_positions[0, 4, :], planet_positions[1, 4, :])
plt.plot(planet_positions[0, 5, :], planet_positions[1, 5, :])
#plt.plot(planet_positions[0, 6, :], planet_positions[1, 6, :])
#plt.legend(["Home", "Target", "Planet 2", "Planet 3", "Planet 4", "Planet 5", "Gas Giant" ])
#plt.legend(["Planet 0(Numerical)", "Planet 5(Numerical)", "Planet 0(Analytical)", "Planet 5(Analytical)"], loc="upper right")
plt.xlabel("AU")
plt.ylabel("AU")
#plt.xlim(-8, 8)
#plt.ylim(-8, 8)
plt.plot(-analytical_orbit(0)[0], -analytical_orbit(0)[1])
plt.plot(-analytical_orbit(5)[0], -analytical_orbit(5)[1])
#plt.plot(analytical_orbit(6)[0], analytical_orbit(6)[1])
plt.legend(["Planet 0(Numerical)", "Planet 5(Numerical)", "Planet 0(Analytical)", "Planet 5(Analytical)"], loc="upper right")
plt.savefig("planet_orbits_analytical.png")
plt.show()