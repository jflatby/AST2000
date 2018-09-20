from __future__ import division
import numpy as np
from ast2000solarsystem_27_v5 import AST2000SolarSystem
import matplotlib.pyplot as plt
import scipy.interpolate as inter
plt.style.use("ggplot")

in_file = open("planet_positions.npy", "rb")
planet_positions, times = np.load(in_file)
print "Positions loaded from file..."

positions_function = inter.interp1d(times, planet_positions)
print "Positions interpolated..."

## Constants
star_system = AST2000SolarSystem(11466)

ms_to_AUyr = 0.000210945021
m_to_AU = 6.68458712226844e-12
s_to_yr = 3.170979198376458591e-8

number_of_planets = star_system.number_of_planets
planet_mass = star_system.mass                                  #[Solar masses]
planet_radius = star_system.radius*1000                         #[m]
semi_major_axes = star_system.a                                 #[AU]

gravitational_constant = 4 * np.pi**2                           #[AU^3 * Yr^(-2) * star masses^(-1)]
star_mass = star_system.star_mass                               #[Solar masses]
star_radius = star_system.star_radius*1000                      #[m]
star_temperature = star_system.temperature                      #[K]
reduced_mass = star_mass*planet_mass/(star_mass + planet_mass)

stefan_boltzmann_const = 5.670373e-8


total_time = 11.4 - 11.3


time_steps_per_year = 70000
total_time_steps = int(time_steps_per_year*total_time)
dt = total_time/total_time_steps
satellite_position = np.zeros((2, total_time_steps))
satellite_velocity = np.zeros((2, total_time_steps))
satellite_acceleration = np.zeros((2, total_time_steps))

##Initial Values
fuel_used_during_launch = 73757.90877855898
launch_duration = 1105.029752 * s_to_yr
launch_time_in_years = 11.3#11.3336428571 #11.27#+ 11.3545
start_time = launch_time_in_years + launch_duration

planet_velocity_at_launch = (positions_function(start_time)[:, 0] - positions_function(launch_time_in_years)[:, 0]) / launch_duration

satellite_mass = 1100 + 100000 - fuel_used_during_launch

def apply_gravity_force(time, sat_pos_rel_star):
    temp_acc = np.array([0.0, 0.0])

    #Star
    r = np.linalg.norm(sat_pos_rel_star)
    temp_acc += -gravitational_constant * star_mass * sat_pos_rel_star / r**3

    #Planets
    for p in range(number_of_planets):
        r_vec = sat_pos_rel_star - positions_function(time)[:, p]
        r = np.linalg.norm(r_vec)
        temp_acc += -gravitational_constant * planet_mass[p] * r_vec / r**3


    return temp_acc

def calculate_launch_time_for_hohmann_transfer():
    angle_wanted = np.pi * (1 - 1 / (2 * np.sqrt(2)) * np.sqrt((semi_major_axes[0] / semi_major_axes[1] + 1) ** 3))
    epsilon = 1e-2
    years = 15

    angle_difference = star_system.omega[1] - star_system.omega[0]
    time_steps_passed = 0
    while abs(angle_wanted-angle_difference) > epsilon and time_steps_passed < time_steps_per_year*years:
        time_steps_passed += 1
        angle_home_planet = np.arctan2(positions_function(time_steps_passed*dt)[1, 0], positions_function(time_steps_passed*dt)[0, 0])
        angle_target_planet = np.arctan2(positions_function(time_steps_passed*dt)[1, 1], positions_function(time_steps_passed*dt)[0, 1])
        angle_difference = angle_target_planet - angle_home_planet

    launch_time_wanted = time_steps_passed*dt
    if abs(angle_wanted-angle_difference) <= epsilon:
        return launch_time_wanted
    else:
        print "couldn't find angle"

def rotate_vector_by(vec, theta):
    rot_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                           [np.sin(theta), np.cos(theta)]])
    return np.dot(rot_matrix, vec)

def compute_planet_velocities(t):
		if t == 0:
			return np.array([star_system.vx0, star_system.vy0])
		elif t == total_time:
			v=(positions_function(start_time + t)-positions_function(start_time + t-dt))/float(dt)
		else:
			v=(positions_function(start_time + t +dt)- positions_function(start_time + t-dt))/float(2*dt)

		return v

print "Pafjlns BOKS: ", compute_planet_velocities(15.8890948351-start_time)[:, 1]

##Function from part 1:
def fuel_left_after_deltav(init_mass, deltav):
    engine_force = 1073357.35327                #[N]
    fuel_loss_per_second = 66.7474270958        #[kg/s]

    exhaust_velocity = engine_force/fuel_loss_per_second
    mf = init_mass * np.exp(-deltav/exhaust_velocity)
    return mf

def orbital_insertion_dv(sat_vel, target_planet, planet_vel, planet_sat_vec):
    sat_vel_rel_to_planet = sat_vel - planet_vel
    planet_sat_distance = np.absolute(np.linalg.norm(planet_sat_vec))

    orbital_vel = np.sqrt(gravitational_constant * planet_mass[target_planet] / planet_sat_distance)

    radial_unit_vector = planet_sat_vec / planet_sat_distance

    orbital_unit_vector = np.array([radial_unit_vector[1], -radial_unit_vector[0]])

    orbital_vel_vector = orbital_vel * orbital_unit_vector

    needed_delta_v = orbital_vel_vector - sat_vel_rel_to_planet

    print "needed deltav: ", needed_delta_v
    return needed_delta_v


"""
Launch time: +11.3
Init vel: + 0.10112
closest approach:  6.92082527303e-05
[  6.51348531e-05   2.33930151e-05]
At time after launch:  4.22702857143
"""
distance_from_planet_center = 10641423.142761754 * m_to_AU
escape_velocity = (10789.272032602184) * ms_to_AUyr  #+  0.0999#0.100680 #0.10112  #11269 m/s

launch_angle_unit_vector = planet_velocity_at_launch/np.linalg.norm(planet_velocity_at_launch)
initial_satellite_position_after_launch_rel_planet = distance_from_planet_center * launch_angle_unit_vector

sat_pos_copy = initial_satellite_position_after_launch_rel_planet.copy()
w = 2.0 * np.pi * star_system.period[0] / 86400.0
theta = w * launch_duration/s_to_yr
initial_satellite_position_after_launch_rel_planet[0] = sat_pos_copy[0] * np.cos(theta) - sat_pos_copy[1] * np.sin(theta)
initial_satellite_position_after_launch_rel_planet[1] = sat_pos_copy[0] * np.sin(theta) + sat_pos_copy[1] * np.cos(theta)

init_rot = np.array((-initial_satellite_position_after_launch_rel_planet[1], initial_satellite_position_after_launch_rel_planet[0]))
final_rot_unit = init_rot / np.linalg.norm(init_rot)
velocity_from_planet_rotation = final_rot_unit * w * distance_from_planet_center

initial_satellite_position_after_launch = positions_function(start_time)[:, 0] + initial_satellite_position_after_launch_rel_planet

initial_satellite_velocity = planet_velocity_at_launch + escape_velocity*launch_angle_unit_vector + velocity_from_planet_rotation
#Overwriter denne med velocity gitt fra mass_needed_launch
initial_satellite_velocity = np.array([1.98825741,  6.58785222]) + 0.09806 * launch_angle_unit_vector
#[1.98825741,  6.58785222
print "Her: ", 0.05635 * launch_angle_unit_vector

satellite_position[:, 0] = initial_satellite_position_after_launch
satellite_velocity[:, 0] = initial_satellite_velocity
satellite_acceleration[:, 0] = apply_gravity_force(start_time, initial_satellite_position_after_launch)

closest_approach = 100
time_of_closest_approach = 0
boost1 = False

def boost(sat_mass, deltav, current_vel):
    new_sat_mass = sat_mass - fuel_needed_for_deltav(sat_mass, deltav)
    new_velocity = current_vel + deltav
    return new_velocity, new_sat_mass


for t in xrange(total_time_steps-1):
    satellite_position[:, t+1] = satellite_position[:, t] + satellite_velocity[:, t] * dt + 0.5 * satellite_acceleration[:, t] * dt ** 2
    satellite_acceleration[:, t+1] = apply_gravity_force(times[t+1]+start_time, satellite_position[:, t+1])
    satellite_velocity[:, t+1] = satellite_velocity[:, t] + 0.5 * (satellite_acceleration[:, t] + satellite_acceleration[:, t+1]) * dt


    satellite_planet_vec = positions_function(times[t+1] + start_time)[:, 1] - satellite_position[:, t+1]
    satellite_planet_distance = np.absolute(np.linalg.norm(satellite_planet_vec))


    if satellite_planet_distance < 0.0001 and boost1 == False:
        print "Performing orbital injection boost at time: ", t*dt, " and position: ", satellite_position[:, t]
        print positions_function(t*dt+start_time)[:, 1]
        deltav = orbital_insertion_dv(satellite_velocity[:, t+1], 1, compute_planet_velocities(t*dt)[:, 1], satellite_planet_vec)
        new_velocity, new_mass = boost(satellite_mass, deltav, satellite_velocity[:, t+1])

        satellite_velocity[:, t+1] = new_velocity
        satellite_mass = new_mass
        boost1 = True


    #if satellite_planet_distance < np.linalg.norm(satellite_position[:, t]) * np.sqrt(planet_mass[1]/(star_mass*10)):
    #    print "Satellite close enough to planet!"


    if satellite_planet_distance < closest_approach:
        closest_approach_vec = satellite_planet_vec
        closest_approach = satellite_planet_distance
        time_of_closest_approach_after_launch = dt*t


    if (dt*t).is_integer() or (dt*t + 0.5).is_integer():
        print (str(dt*t) + " years calculated")

print "closest approach: ", closest_approach
print closest_approach_vec
print "At time after launch: ", time_of_closest_approach_after_launch
print "distance from surface: ", closest_approach - planet_radius[1]*m_to_AU

## Plot
print "Starting plotting..."
plt.plot(satellite_position[0, :], satellite_position[1, :], label="Satellite")

#plt.scatter(5.58197585547 ,  -2.87298510355 , s=10)

#sat_planet_vec = np.array([5.58197585547 ,  -2.87298510355 ])-positions_function(20)[:, 1]
#sat_planet_unit = sat_planet_vec/np.linalg.norm(sat_planet_vec)
#print np.linalg.norm(sat_planet_vec)

#orbital_insertion_dv(np.array([1.73894168678 , -3.8306764623]), 1, compute_planet_velocities(total_time)[:, 1], sat_planet_vec)

for p in xrange(2):

    plt.plot(positions_function(times[int(start_time/dt):int(start_time/dt)+total_time_steps])[0, p], positions_function(times[int(start_time/dt):int(start_time/dt)+total_time_steps])[1, p], label=("Planet ", p))

plt.legend(["Satellite", "Home planet", "Target Planet"])
plt.xlabel("AU")
plt.ylabel("AU")
plt.savefig("hohmann_transfer_sim.png")
plt.show()

init_sat_pos = positions_function(launch_time_in_years)[:, 0] + launch_angle_unit_vector*planet_radius[0]*m_to_AU

print init_sat_pos, launch_time_in_years, satellite_position[:, 0]

#star_system.engine_settings(3.4021029823078977e-09, 3.15498195925e14, 63200999999999.992, 100000, 1105.029752, init_sat_pos, launch_time_in_years)
#print star_system.mass_needed_launch(satellite_position[:, 0], test=True)


print initial_satellite_velocity
# SKal f array([ 1.90295741,  6.61099202])
print "mass after launch: ", satellite_mass
satellite_mass = fuel_needed_for_deltav(satellite_mass, np.linalg.norm([0.01895096, 0.05306773])/ms_to_AUyr)
print "mass after first boost: ", satellite_mass
satellite_mass = fuel_needed_for_deltav(satellite_mass, np.linalg.norm([-0.46930108, -0.02598787])/ms_to_AUyr)
print "mass after orbital insertion: ", satellite_mass
