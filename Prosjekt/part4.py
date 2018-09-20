from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from ast2000solarsystem_27_v4 import AST2000SolarSystem
from PIL import Image
import scipy.interpolate as inter
plt.style.use("ggplot")
star_system = AST2000SolarSystem(11466)

ms_to_AUyr = 0.000210945021

in_file = open("planet_positions.npy", "rb")
planet_positions, times = np.load(in_file)
print "Positions loaded from file..."

positions_function = inter.interp1d(times, planet_positions)
print "Positions interpolated..."

with open('himmelkule.npy') as infile:
    himmelkulen = np.load(infile)


def compute_reference_sphere():
    x_pixels = 640
    y_pixels = 480
    field_of_view = 1.22173

    max_axis_value = 2 * np.sin(field_of_view / 2) / (1 + np.cos(field_of_view / 2))
    picture = np.zeros((y_pixels, x_pixels, 3), dtype=np.uint8)
    x = np.linspace(-max_axis_value, max_axis_value, x_pixels)
    y = np.linspace(-max_axis_value, max_axis_value, y_pixels)

    theta_0=np.pi/2.0
    pictures=np.zeros(shape=(360, int(y.shape[0]), int(x.shape[0]), 3), dtype=np.uint8)
    xx, yy=np.meshgrid(x,y)
    rho=np.sqrt(xx**2+yy**2)
    c=2.0*np.arctan(rho/2.0)
    theta=theta_0-np.arcsin(yy*np.sin(c)/rho) #Is sin(theta_0) not simply zero??

    for phi_0 in range(0, 360):
        phi_in_rad=np.deg2rad(phi_0)
        phi=phi_in_rad+np.arctan(xx*np.sin(c)/(rho*np.cos(c)))
        for k in range(len(x)):
            for j in range(len(y)):
                pixnum=AST2000SolarSystem.ang2pix(theta[j,k], phi[j,k])
                temp=himmelkulen[pixnum]
                pictures[phi_0, j, k, :]=[temp[2], temp[3], temp[4]]
        print "Done with phi: ", phi_0
    np.save("Reference_sphere.npy", pictures)

def generate_projection(phi_midpoint):
    theta_midpoint = np.pi / 2
    x_pixels = 640
    y_pixels = 480
    field_of_view = 1.22173

    max_axis_value = 2 * np.sin(field_of_view / 2) / (1 + np.cos(field_of_view / 2))
    picture = np.zeros((y_pixels, x_pixels, 3), dtype=np.uint8)
    x_values = np.linspace(-max_axis_value, max_axis_value, x_pixels)
    y_values = np.linspace(-max_axis_value, max_axis_value, y_pixels)


    for xpix in xrange(x_pixels):
        for ypix in xrange(y_pixels):
            rho = np.linalg.norm([x_values[xpix], y_values[ypix]])
            c = 2 * np.arctan2(rho, 2)


            theta = np.pi / 2 - np.arcsin(np.cos(c) * np.cos(theta_midpoint) + (y_values[ypix] * np.sin(c) * np.sin(theta_midpoint)) / rho)
            phi = phi_midpoint + np.arctan2(x_values[xpix] * np.sin(c), rho * np.sin(theta_midpoint) * np.cos(c) - y_values[ypix] * np.cos(theta_midpoint) * np.sin(c))


            pixel_number = AST2000SolarSystem.ang2pix(theta, phi)
            rgb = [himmelkulen[pixel_number][2], himmelkulen[pixel_number][3], himmelkulen[pixel_number][4]]

            picture[ypix, xpix, :] = rgb

    return picture


def use_projection():
    x_pixels = 640
    y_pixels = 480
    phi = 0
    degrees = 360

    complete_picture = np.zeros((degrees, y_pixels, x_pixels, 3), dtype=np.uint8)

    for degree in xrange(degrees):
        complete_picture[degree, :, :, :] = generate_projection(degree)
        print degree

    np.save('sky.npy', complete_picture)



def find_best_orientation_angle(picture):
    with open('sky.npy', 'rb') as infile:
        sky = np.load(infile)

    degrees = 360
    least_error = 100000000
    least_error_deg = 0

    for degree in xrange(degrees):


        if degree % 10 == 0:
            test = Image.fromarray(sky[degree, :, :, :])
            test.save("test%i.png" %(degree))


        pic_from_file = sky[degree, :, :, :]

        #error = np.linalg.norm(pic_from_file - picture)
        diff_img = (picture-pic_from_file)**2
        error = np.sum(np.sum(np.sum(diff_img)))

        if error < least_error:
            least_error = error
            least_error_deg = degree


    return least_error_deg

def get_radial_velocities(measured_shift1, measured_shift2):
    """
    Reference star 1, at phi = 96.394216 degrees,
    has delta lambda of  -0.002651878563 nm in the home star rest frame.
    Reference star 2, at phi = 56.252808 degrees,
    has delta lambda of   0.008295625752 nm in the home star rest frame.
    """


    speed_of_light = 63197.7909261                                  # [AU/Yr]
    h_alpha_spectral_line = 656.3*1e-9                                   # [nanometer] rest wavelength

    phi_star_1 = 96.394216
    phi_star_2 = 56.252808

    radial_velocity_1 = measured_shift1 * speed_of_light / h_alpha_spectral_line
    radial_velocity_2 = measured_shift2 * speed_of_light / h_alpha_spectral_line

    return radial_velocity_1, radial_velocity_2


def get_satellite_velocity():
    """
    The radial velocity of the reference star relative to the satellite is v_rel = v_refstar - v_sat.
    Where v_sat is the velocity of the satellite in the direction of the reference star.
    v_sat in the direction phi_1 will be called v_1, and in direction phi_2
    """
    placeholder_lambda1, placeholder_lambda2 = -0.055629809391*1e-9, -0.040744894307*1e-9

    phi_star_1 = np.deg2rad(96.394216)
    phi_star_2 = np.deg2rad(56.252808)

    star_1_radial, star_2_radial = get_radial_velocities(-0.002651878563*1e-9, 0.008295625752*1e-9)
    relative_velocity_sat_star_1, relative_velocity_sat_star_2 = get_radial_velocities(placeholder_lambda1, placeholder_lambda2)
    relative_velocity_sat = np.array([star_1_radial - relative_velocity_sat_star_1, star_2_radial - relative_velocity_sat_star_2])

    factor = 1 / np.sin(phi_star_2 - phi_star_1)

    inverse_of_radial_matrix = np.array(([np.sin(phi_star_2), -np.sin(phi_star_1)], [-np.cos(phi_star_2), np.cos(phi_star_1)]))

    satellite_velocity = factor * (np.dot(inverse_of_radial_matrix, relative_velocity_sat))

    return satellite_velocity


number_of_planets = 7
distances_to_planets = np.load("pos.npy")


def get_satellite_position(distances, current_time):
    """
    Calculates the satellite position from a list of distances the current time.
    pretty much looks at the circles with radius[i] and checks for intersections
    """

    radius_star = distances[-1]
    radii_planets = distances[:-1]
    planet_positions = positions_function(current_time).T

    angle_array = np.linspace(0, 2*np.pi, 10000)
    star_distance_circle = radius_star * np.array((np.sin(angle_array), np.cos(angle_array)))
    plt.plot(star_distance_circle[0], star_distance_circle[1], color='k')


    for planet in xrange(number_of_planets):
        x0 = planet_positions[planet, 0]
        y0 = planet_positions[planet, 1]
        radius = radii_planets[planet]

        A = (x0**2 + y0**2 + radius_star**2 - radius**2) / (2 * x0)
        B = y0 / x0
        C = np.sqrt(radius_star**2 * (B**2 + 1) - A**2)

        y = [(A*B + C) / (B**2 + 1), (A*B - C) / (B**2 + 1)]
        x = [A - i * B for i in y]

        for i in xrange(2):
            plt.scatter(x[i], y[i])

        circle = np.array((radius * np.sin(angle_array) + x0, radius * np.cos(angle_array) + y0))
        plt.plot(circle[0], circle[1])

    plt.scatter(0,0,c='y', s=100)
    plt.axis('equal')
    plt.savefig("triangulation.png")
    plt.show()

#use_projection()
#compute_reference_sphere()

img = Image.open("find_orient.png")
pixels = np.array(img)
img.close()

#print find_best_orientation_angle(pixels)
get_satellite_position(distances_to_planets, 11.5)