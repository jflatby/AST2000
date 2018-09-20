
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.patches as patches
import mpl_toolkits.mplot3d.art3d as art3d
import numpy as np
from ast2000solarsystem_27_v3 import AST2000SolarSystem

np.random.seed(1)
solar_system_seed = 11466


H2_molar_mass = 2.01588            #[g/mol]
avogadros_const = 6.02214179e23
boltzmann_const = 1.380648e-23
particle_mass = H2_molar_mass/(avogadros_const * 1e3) #[kg]
gravitational_constant = 6.67408e-11  # [m**3 * kg**-1 * s**-2]

satellite_mass = 1100  # [kg]

solar_system = AST2000SolarSystem(solar_system_seed)
home_planet_mass = 1.98855e30 * solar_system.mass[0]  # [kg]
home_planet_radius = solar_system.radius[0] * 1000  # [m]
home_planet_gravity = gravitational_constant * home_planet_mass / home_planet_radius ** 2
home_planet_escape_velocity = np.sqrt(2*gravitational_constant*home_planet_mass/home_planet_radius)

def simulate_box(number_of_particles, Nt, dt, box_dim, T):
    """
    Calculates the motion of the particles used for fuel calculations
    variables:
    n  =        number of particles
    Nt =        Number of time steps

    dt =        Time step
    box_dim =   Dimension of the box
    """

    hole_size = box_dim / 2
    hole = (box_dim / 2 - hole_size / 2, box_dim / 2 + hole_size / 2)

    sigma = np.sqrt(boltzmann_const * T / particle_mass)     # Standard deviation for initial gaussian velocity distribution
    mean = 0  # Mean for initial gaussian distribution (no bulk flow)

    positions = np.zeros((number_of_particles, 3, Nt))  # Position vector

    velocities = np.zeros((number_of_particles, 3))  # Velocity vector

    """ Setting initial positions and velocities. Positions are selected from a uniform distribution while
    velocities are selected from a gaussian distribution with a mean and deviation described in part1. """

    for j in range(number_of_particles):
        for k in range(3):
            positions[j, k, 0] =  np.random.uniform(0, box_dim)
            velocities[j, k] =  np.random.normal(mean, sigma)

    partition_momentum = 0.0  # Momentum value scalar in x-direction
    number_of_escapes = 0  # Counter for number of particles escaped

    for i in range(Nt - 1):  # Time loop
        positions[:, :, i+1] = positions[:, :, i] + velocities[:, :] * dt  # Move ALL particles one step

        outside_box_mask1 = positions[:, :, i + 1] < 0
        outside_box_mask2 = positions[:, :, i + 1] > box_dim
        velocities[outside_box_mask1] *= -1
        velocities[outside_box_mask2] *= -1

        x_mask = positions[:, 0, i + 1] < 0
        y_mask = np.logical_and(positions[:, 1, i + 1] > hole[0], positions[:, 1, i + 1] < hole[1])
        z_mask = np.logical_and(positions[:, 2, i + 1] > hole[0], positions[:, 2, i + 1] < hole[1])

        xy_mask = np.logical_and(x_mask, y_mask)
        xyz_mask = np.logical_and(xy_mask, z_mask)

        particles_in_hole = np.sum(xyz_mask)
        number_of_escapes += particles_in_hole

        partition_momentum += particle_mass * particles_in_hole * sum(velocities[xyz_mask, 0])

    return partition_momentum, number_of_escapes

def launch_sim(no_boxes, v_end, A, no_particles_sec):
    """
    Simulates the launch.

    Variables:
    no_boxes:           the number of fuel boxes you need
    v_end:              the desired velocity (relative to surface) at the end of launch
    A:                  dp/dt for the box described in first method
    no_particles_sec:   number of particles pr second leaving the box in first method

    returns:    How much time it takes to achieve the boost and how much fuel you have used up.
    """

    A = A * no_boxes  # (dp/dt) for whole engine
    B = 0.  # (insert value)  # dm/dt for the box described in first method
    v = 0.  # initial velocity relative to surface of planet
    time = 0.  # initialize time variable

    T = 20. * 60  # Total time, 20 minutes
    Nt = 10000  # Number of time steps
    dt = float(T) / Nt  # Time step

    initial_fuel_mass = 100  # Calculate/set fuel mass
    M = self.mass_sat + initial_fuel_mass  # Total mass

    for i in xrange(Nt):
        """ Here you need to update the new value of the velocity aswell as the new total mass value M """
        # v += boost*dt - gravity ...
        # M -= ...
        time += dt
        if M < self.mass_sat:
            """ Some sort of error message telling you you're out of fuel """
        elif v < 0:
            """ Gravity stronger than engine """
        elif v >= v_end:
            """ Boost is successful. Save values and end method """
            # fuel_needed = ...
            print "Boost Succesful"
            return fuel_needed, time

    return 0., 0.  # returns 0 because the boost was not successful.

def plot(number_of_particles, Nt, dt, box_dim, T):
    """ A method for plotting and animating a simulation of particles in box """

    def update_lines(num, dataLines, lines):
        for line, data in zip(lines, dataLines):
            line.set_data(data[0:2, num - 1:num])
            line.set_3d_properties(data[2, num - 1:num])
        return lines

    # Attach 3D axis to the figure
    fig = plt.figure()
    ax = p3.Axes3D(fig)

    frames = 100

    # Run the actual simulation
    x_momentum, datax = simulate_box(number_of_particles, Nt, dt, box_dim, T)
    lines = []
    data = []
    for i in range(number_of_particles):
        data.append([datax[i]])  # wrap data inside another layer of [], needed for animation!
        lines.append([ax.plot(data[i][0][0, 0:1], data[i][0][1, 0:1], data[i][0][2, 0:1], 'o')[0]])

    # Set the axes properties
    ax.set_xlim3d([0.0, box_dim])
    ax.set_xlabel('X')

    ax.set_ylim3d([0.0, box_dim])
    ax.set_ylabel('Y')

    ax.set_zlim3d([0.0, box_dim])
    ax.set_zlabel('Z')

    ax.set_title('Particle Animation')

    hole = patches.Circle((box_dim / 2, box_dim / 2), box_dim / 4)
    hole.set_facecolor('black')
    ax.add_patch(hole)
    art3d.pathpatch_2d_to_3d(hole, z=0, zdir="z")

    ani = [i for i in xrange(number_of_particles)]
    for i in ani:
        ani[i] = animation.FuncAnimation(fig, update_lines, frames, fargs=(data[i], lines[i]),
                                         interval=50, blit=False)

    plt.show()


#print simulate_box(box_dim=1e-6, number_of_particles=100000, Nt=1000, dt=1e-12, T=10000)
plot(100, )

