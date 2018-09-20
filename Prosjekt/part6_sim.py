from __future__ import division
import numpy as np
from ast2000solarsystem_27_v5 import AST2000SolarSystem
import matplotlib.pyplot as plt
import scipy.interpolate as inter
from scipy.constants import gravitational_constant, Boltzmann, m_p
plt.style.use("ggplot")
in_file = open("planet_positions.npy", "rb")
planet_positions, times = np.load(in_file)
print "Positions loaded from file..."

positions_function = inter.interp1d(times, planet_positions)
print "Positions interpolated..."

rho_values, r_values = np.load("rho_values.npy")
print r_values[0], r_values[-1]
rho_func = inter.interp1d(r_values, rho_values)

## Constants
star_system = AST2000SolarSystem(11466)
ms_to_AUyr = 0.000210945021
m_to_AU = 6.68458712226844e-12
s_to_yr = 3.170979198376458591e-8
solar_mass_to_kg = 1.99e30

number_of_planets = star_system.number_of_planets
planet_mass = star_system.mass                                  #[Solar masses]
planet_radius = star_system.radius*1000                         #[m]
semi_major_axes = star_system.a                                 #[AU]
star_mass = star_system.star_mass                               #[Solar masses]
star_radius = star_system.star_radius*1000                      #[m]
star_temperature = star_system.temperature                      #[K]
reduced_mass = star_mass*planet_mass/(star_mass + planet_mass)
satellite_mass = 23412.1651052
gamma = 1.4
surface_temp = 271
rho0 = star_system.rho0[1]
lander_mass = 90
gravity = gravitational_constant * star_system.mass[1] * solar_mass_to_kg / (star_system.radius[1] * 1000)**2
mu2 = 17.02 # mean molecular weight
scale_height = Boltzmann * surface_temp / (2 * m_p * mu2 * gravity)

lander_area = 0.3
parachute_area = 22.2

planet_rotational_speed = (2 * np.pi * planet_radius[1])/(star_system.period[1]*24*60*60)

mu = gravitational_constant*planet_mass[1]*solar_mass_to_kg

stefan_boltzmann_const = 5.670373e-8

total_time = 20000#11300 #0.2 * 86400 #[seconds]

time_steps_per_second = 10
total_time_steps = int(time_steps_per_second*total_time)
dt = total_time/total_time_steps

satellite_position = np.zeros((2, total_time_steps))
satellite_velocity = np.zeros((2, total_time_steps))
satellite_acceleration = np.zeros((2, total_time_steps))

##Initial Values
start_time = 15.88905 #[Years]
def compute_planet_velocities(t): #Time after start time
    dt2 = 1.42877553936e-05
    if t == 0:
        return np.array([star_system.vx0, star_system.vy0])
    elif t == total_time:
        v=(positions_function(t)-positions_function(t-dt2))/float(dt2)
    else:
        v=(positions_function(t +dt2)- positions_function(t-dt2))/float(2*dt2)

    return v

init_planet_vel_rel_star = compute_planet_velocities(start_time)[:, 1] / ms_to_AUyr
init_planet_pos_rel_star = positions_function(start_time)[:, 1] / m_to_AU
init_sat_vel_rel_star = np.array([1.26964060678, -3.8566643323]) / ms_to_AUyr
init_sat_pos_rel_star = np.array([-6.24349802127, -0.273990555621]) / m_to_AU
init_sat_pos_rel_planet = init_sat_pos_rel_star - init_planet_pos_rel_star
init_sat_vel_rel_planet = init_sat_vel_rel_star - init_planet_vel_rel_star


## Nu er planeten i origo

def apply_gravity_force(sat_pos):
    r_vec = sat_pos
    r = np.linalg.norm(r_vec)
    acc = -gravitational_constant * (planet_mass[1]*solar_mass_to_kg) * r_vec / r**3
    return acc

def find_r_isothermal():
    r_inv = (1.0 / star_system.radius[1] * 1000) - (surface_temp / 2) * (gamma / (gamma - 1)) * (Boltzmann / (gravitational_constant * planet_mass[1]*solar_mass_to_kg * mu2 * m_p))
    return r_inv**(-1)

parachute = False

def apply_drag_force(sat_pos, sat_vel):
    r = np.linalg.norm(sat_pos)
    if r-planet_radius[1] > 100000:
        return np.array([0, 0])

    radial_unit_vec = sat_pos / r
    angular_unit_vec = np.array([-radial_unit_vec[1], radial_unit_vec[0]])
    #print "sat_vel: ", sat_vel - planet_rotational_speed*angular_unit_vec
    sat_vel_rel_atmosphere = sat_vel #- planet_rotational_speed*angular_unit_vec
    #print "planet rot: ", planet_rotational_speed

    if r >= 6128701.23927 or r <= planet_radius[1]:
        rho = 0
    else:
        rho = rho_func(r)


    if r-planet_radius[1] > 5000:
        drag_force = (0.5 * rho * np.linalg.norm(sat_vel_rel_atmosphere)**2 * lander_area)
    else:
        #print "parachute"
        drag_force = (0.5 * rho * np.linalg.norm(sat_vel_rel_atmosphere)**2 * parachute_area)

    velocity_unit_vec = sat_vel_rel_atmosphere / np.linalg.norm(sat_vel_rel_atmosphere)
    acc = drag_force / lander_mass * -velocity_unit_vec

    if np.isnan(acc).all():
        acc=np.array([0, 0])

    return acc

satellite_position[:, 0] = init_sat_pos_rel_planet
satellite_velocity[:, 0] = init_sat_vel_rel_planet
satellite_acceleration[:, 0] = apply_gravity_force(init_sat_pos_rel_planet)

boost3_time = 0 # Bruke for aa velge hvor vi lander?

deltav_1 = 0
boost2_time = -1
r1 = 0
r2 = 0

for t in xrange(total_time_steps-1):
    satellite_position[:, t+1] = satellite_position[:, t] + satellite_velocity[:, t] * dt + 0.5 * satellite_acceleration[:, t] * dt ** 2
    F_d = apply_drag_force(satellite_position[:, t], satellite_velocity[:, t])
    #print F_d
    satellite_acceleration[:, t+1] = apply_gravity_force(satellite_position[:, t+1]) + F_d
    satellite_velocity[:, t+1] = satellite_velocity[:, t] + 0.5 * (satellite_acceleration[:, t] + satellite_acceleration[:, t+1]) * dt


    if np.linalg.norm(satellite_position[:, t]) <= planet_radius[1]:
        #print "Crash"
        satellite_position[:, t+1] = satellite_position[:, t]
        satellite_velocity[:, t+1] = satellite_velocity[:, t]
        satellite_acceleration[:, t+1] = satellite_acceleration[:, t]

    if (t*dt)%(total_time/5) == 0:
        print t*dt

    if t == int(2000/dt):
        #plt.scatter(satellite_position[0, t], satellite_position[1, t])
        #dot1 = satellite_position[:, t]

        r1 = np.linalg.norm(satellite_position[:, t+1])
        r2 = 6150000#6030000
        print r2 - planet_radius[1]
        deltav_1 = np.sqrt(mu / r1) * (np.sqrt((2 * r2) / (r1 + r2)) - 1)
        radial_unit_vec = satellite_position[:, t+1]/np.linalg.norm(satellite_position[:, t+1])
        angular_unit_vec = -np.array([-radial_unit_vec[1], radial_unit_vec[0]])
        satellite_velocity[:, t+1] += deltav_1*angular_unit_vec
        boost2_time = int(2000 + (np.pi * np.sqrt((r1+r2)**3/(8*mu))))
        print "boost1_time: ", 15.88905 + 2000 * s_to_yr
        print "boost2_time: ", 15.88905 + boost2_time * s_to_yr
        print "deltav1: ", deltav_1 * angular_unit_vec * ms_to_AUyr

    if t == int(boost2_time/dt):
        #plt.scatter(satellite_position[0, t], satellite_position[1, t])
        #dot2 = satellite_position[:, t]

        deltav_2 = np.sqrt(mu / r2) * (1 - np.sqrt((2 * r1) / (r1 + r2)))
        radial_unit_vec = satellite_position[:, t+1] / np.linalg.norm(satellite_position[:, t+1])
        angular_unit_vec = -np.array([-radial_unit_vec[1], radial_unit_vec[0]])
        satellite_velocity[:, t+1] += deltav_2*angular_unit_vec
        print "deltav2: ", deltav_2*angular_unit_vec * ms_to_AUyr

    boost3_time_after2 = 6000
    boost3_deltav = 50
    if t == int(boost2_time/dt + boost3_time_after2/dt):
        boost3_time = int(boost2_time/dt + boost3_time/dt)
        vel_unit_vec = satellite_velocity[:, t]/np.linalg.norm(satellite_velocity[:, t])
        satellite_velocity[:, t+1] += (boost3_deltav)*-vel_unit_vec
        print "deltav3: ", (boost3_deltav)*-vel_unit_vec


    if np.linalg.norm(satellite_position[:, t]) - planet_radius[1] <= 100000 and parachute == False:
        print "parachute time: ", t*dt - boost2_time
        parachute = True

    if t == int(((15.89-start_time)/s_to_yr)/dt):
        plt.scatter(satellite_position[0, t], satellite_position[1, t])


#plt.plot([dot1[0], dot2[0]], [dot1[1], dot2[1]])

print "period at 7000000: ", 2*np.pi * np.sqrt(r2**3 / mu)

print "Starting plotting..."
plt.plot(satellite_position[0, :], satellite_position[1, :], label="Satellite")
print boost3_time*dt
#Planet
fig = plt.gcf()
ax = fig.gca()
circle = plt.Circle((0.0, 0.0), planet_radius[1])
circle2 = plt.Circle((0.0, 0.0), 200000, color="black")

ax.add_artist(circle)
ax.add_artist(circle2)


plt.xlim(-15000000, 15000000)
plt.ylim(-15000000, 15000000)

plt.legend(["Satellite"])
plt.xlabel("m")
plt.ylabel("m")
plt.show()


"""
6000000m
deltav1:  [-0.16080915  0.133695  ]
deltav2:  [ 0.19414536 -0.16316268]

"""
