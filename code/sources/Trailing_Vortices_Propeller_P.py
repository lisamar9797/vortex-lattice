"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine calculates the induced velocities (coefficient)
from the transition wake and from the semi-infinite helicoidal vortex at a
single point for all propeller blades. Specifically, it addresses the selected
trailing vortex (n_tral_vortex) - options 1, 2, 3, 4 with Msp set to 3.
"""


import numpy as np
import sources.Variables as Var
from sources.De_Jong_P import De_Jong
from sources.Biot_Savart_Propeller_P import Biot_Savart_Propeller
from sources.Helix_P import Helix


def Trailing_Vortices_Propeller(n_tral_vortex, px, py, pz):
    Points_Trans_Wake_P = np.loadtxt("output/Propeller_Points_Trans_Wake.txt", skiprows= 1, usecols= (1,2,3))
    Grid_Points_P = np.loadtxt("output/Propeller_Grid_Points.txt")
    N_Bound_Vortex_P = np.loadtxt("output/Propeller_N_Bound_Vortex.txt",dtype= 'int')
    N_Bound_Vortex_P = N_Bound_Vortex_P.reshape((Var.Msp+1, 1))

    U_x = 0.0         # Initialization of the variable U_x
    U_y = 0.0         # Initialization of the variable U_y
    U_z = 0.0         # Initialization of the variable U_z

    # INDUCED VELOCITIES FROM SEMI-INFINITE HELICOIDAL VORTEX

    pyy = py
    pzz = pz

    k_2 = (n_tral_vortex + n_tral_vortex * Var.N_P_L)+Var.N_P_L
    # k_2 is the location in Points_Trans_Wake_P for the last point of that trailing vortex
    k_0 = k_2 - Var.N_P_L
    # k_0 is the location in Points_Trans_Wake_P for the first point of that trailing vortex (T.E.)
    k_2 = int(k_2)
    k_0 = int(k_0)
    n_tral_vortex = int (n_tral_vortex)

    x_1 = - Points_Trans_Wake_P[k_2,0]
    # First point for the semi-infinite helicoidal vortex (x) The - is because in infv the vortex starts at -infinity and stops at -x

    r_1 = Points_Trans_Wake_P[k_2,1]    # First point for the semi-infinite helicoidal vortex (Radius)

    p_1 = Points_Trans_Wake_P[k_2,2]    # First point for the semi-infinite helicoidal vortex (pitch)

    x_T_E =	Points_Trans_Wake_P[k_0,0]
    r_T_E = Points_Trans_Wake_P[k_0,1]
    p_T_E = Points_Trans_Wake_P[k_0,2]

    delta_theta = 2*np.pi/float(Var.Z_Blade_P)

    for i in range(Var.Z_Blade_P):

        y_T_E = Grid_Points_P[N_Bound_Vortex_P[n_tral_vortex,0],1]
        theta = np.arcsin(- y_T_E / r_T_E) # Theta for the T.E. point
        xki = theta - 2*np.pi* x_T_E/p_T_E # Phase angle phi

        U_xx, U_yy, U_zz = De_Jong (x_1,r_1,p_1,xki,px, pyy, pzz)
        # This subroutine calculates the induced velocity for the ultimate wake that starts in x1,r1,z1 (De Jong) for the point px,pyy,pzz

        U_x = U_x + U_xx
        U_y = U_y + U_yy
        U_z = U_z + U_zz

        theta_blade = (float(i+1)*delta_theta)# This is the angle of the blade that induces velocity on the reference blade

        if (theta_blade < np.pi):
            pyy = py*np.cos(theta_blade) + pz*np.sin(theta_blade)
            pzz = pz*np.cos(theta_blade) - py*np.sin(theta_blade)
            # In order to calculate the induced velocity in the point px,py,pz from the semi-infinite
            # helicoidal vortices we do not change the helix (which is always located on the reference blade),
            # but we change the location of the point and we keep costant the relative distance between the semi-infinite helicoidal vortex
            # and the point (the rotation depends on the location of the blade that induces velocity)
        else:
            pyy = py*np.cos(2*np.pi-theta_blade) - pz*np.sin(2*np.pi-theta_blade)
            pzz = pz*np.cos(2*np.pi-theta_blade) + py*np.sin(2*np.pi-theta_blade)

    # INDUCED VELOCITIES FROM THE TRANSITION WAKE

    d_r = 0.0
    d_p = 0.0

    x_1 = Points_Trans_Wake_P[k_2,0]

    for i in range(Var.N_P_L):
            k2_i = k_2 - (i+1)

            x_2 = Points_Trans_Wake_P[k2_i,0]   # Second point for the selected side of the transition wake (x)
            r_2 = Points_Trans_Wake_P[k2_i,1]	# Second point for the selected side of the transition wake (radius)
            p_2 = Points_Trans_Wake_P[k2_i,2] 	# Second point for the selected side of the transition wake (pitch)

            a_r,b_r,c_r = Helix(x_1,x_2,r_1,r_2,d_r)  # Calculates the coefficients a,b,c used in the polynomium
                                                      # for the radius for that side of the transition wake

            a_p,b_p,c_p = Helix(x_1,x_2,p_1,p_2,d_p)  # Calculates the coefficients a,b,c used in the polynomium
                                                      # for the pitch for that side of the transition wake

            delta_x = (x_2-x_1) / float(Var.sub_interv)

            x_11 = x_1	 # First value of the element line (x) of the transition wake
            r_11 = r_1	 # First value of the element line (radius) of the transition wake
            p_11 = p_1	 # First value of the element line (pitch) of the transition wake

            theta_1 = (xki + (2*np.pi) * x_11) / p_11 # First value of the element line (theta) of the transition wake

            y_11 = - r_11 * np.sin(theta_1)	# First value of the element line (y) of the transition wake
        	# Theta is positive in the other direction (sin(-theta) = - sin(theta))

            z_11 = r_11 * np.cos(theta_1)   # First value of the element line (z) of the transition wake


            for j in range (Var.sub_interv):

                x_12 = x_11 + delta_x	                # Second value of the element line (x) of the transition wake
                r_12 = a_r*(x_12**2) + b_r*x_12 + c_r	# Second value of the element line (radius) of the transition wake
                p_12 = a_p*(x_12**2) + b_p*x_12 + c_p	# Second value of the element line (pitch) of the transition wake

                theta_2 = xki + (2*np.pi) * x_12 / p_12 # Second value of the element line (theta) of the transition wake
                y_12 = - r_12 * np.sin(theta_2)			# Second value of the element line (y) of the transition wake
                z_12 = r_12 * np.cos(theta_2)			# Second value of the element line (z) of the transition wake

                U_x_w,U_y_w,U_z_w = Biot_Savart_Propeller (Var.Z_Blade_P,x_11,y_11,z_11,x_12,y_12,z_12,px,py,pz)

                U_x = U_x + U_x_w
                U_y = U_y + U_y_w
                U_z = U_z + U_z_w

                x_11 = x_12
                y_11 = y_12
                z_11 = z_12

            x_1 = x_2
            r_1 = r_2
            p_1 = p_2

            d_r = 2 * a_r * x_2 + b_r
            d_p = 2 * a_p * x_2 + b_p

    return (U_x, U_y, U_z)

