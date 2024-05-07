'''
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This module contains the fixed variables.
'''


# Density of water (kg/m^3)
rho = 1025
# Velocity of the ship (m/s)
V_Ship = 10
# Convergence criteria
epsi = 0.0001
# Number of subdivisions of the input values
N_Iter = 500
# Total required thrust (N)
Tr = 959477.23


# Preserved total required thrust for calculations
Tr_P = Tr
# Number of panels (spanwise)
Msp = 5
# Number of panels (chordwise)
Nch = 5
# Flat plate coefficient (0: pure rooftop, 0.5: half rooftop, 1: pure flat
# plate)
cny = 0.5

# Angular velocity (rad/s) - Propeller
Omega = 11.780
# Radius (m) - Propeller
Rad_P = 3
# Radius for the hub - Propeller
R_Hub_P = 0.2 * Rad_P
# Number of blades - Propeller
Z_Blade_P = 5
# Skin friction drag coefficient
Skin_Coeff = 0.000
# Number of straight line vortices (Transition Wake)
N_P_L = 5
# Number of subintervals for each line of the transition wake
sub_interv = 60
