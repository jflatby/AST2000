from ast2000solarsystem_27_v5 import AST2000SolarSystem
import numpy as np

star_system = AST2000SolarSystem(11466)

star_system.engine_settings(3.4021029823078977e-09, 3.15498195925e14, 63200999999999.992, 100000, 1120,  [4.18992917, -1.46374675], 11.3)
star_system.mass_needed_launch([ 4.18998484, -1.46357376], test=True)

star_system.send_satellite("instructions.txt")
star_system.land_on_planet(1, "land_on_planet.txt")
#53
#[ 0.69970958  5.21180747]
#4.42711 -0.435386


#Forventet pos [-6.24535362 -0.24178915] at 4.58078571429

#boost 11.300040 0.03297837 0.0923482


#Orbital injection
"""
Recorded current satellite position and velocity (at time  15.895 ):
Position = ( -6.24204829104 ,  -0.297574728556 )
Velocity = ( -0.782630279242 ,  -4.06375394709 )
"""

#planet_vel at orbital injection = [ 0.21120866 -3.93634163]
#planet positions
# ##[-6.24549357 -0.23863241]

#boost 15.889160 -0.03804509 -0.0022238
#video 15.9 1


"""
launch
boost 11.300036   0.01895096  0.05306773
video 15.887 1
boost 15.88905 -0.46930108 -0.02598787
video 15.8893 1
boost 15.88906  -0.26733109 -0.01562596
boost 15.8890948351 0.26733109  0.01562596
"""
