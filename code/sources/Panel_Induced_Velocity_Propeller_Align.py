"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine calculates the induced velocities (coefficient)
from a panel in the point (px,py,pz) for all the blades of the propeller
without including the bound vortex.

Parameters:
 - n_pnl   : number of the panel that induces velocity in the point
             (px,py,pz)
 - mpnl    : number of the panel that containes the point px,py,pz.
 - msid    : number of the side that containes the point px,py,pz.
"""


import numpy as np
import sources.Variables as Var
from sources.Biot_Savart_Propeller_P import Biot_Savart_Propeller


def Panel_Induced_Velocity_Propeller_Align(n_pnl, msid, mpnl, px, py, pz):

    Grid_Points_P = np.loadtxt("output/Propeller_Grid_Points.txt")
    N_Panel_P = np.loadtxt("output/Propeller_Numeration_Panel.txt", dtype='int')

    # DECLARATION OF VARIABLES

    U_x = 0     # Initialization of the variable U_x
    U_y = 0     # Initialization of the variable U_y
    U_z = 0     # Initialization of the variable U_z

    delta_theta = 2*np.pi/float(Var.Z_Blade_P)

    x_10 = Grid_Points_P[N_Panel_P[n_pnl,0],0] # X value for the first point of the chosen panel of the propeller
    y_10 = Grid_Points_P[N_Panel_P[n_pnl,0],1] # Y value for the first point of the chosen panel of the propeller
    z_10 = Grid_Points_P[N_Panel_P[n_pnl,0],2] # Z value for the first point of the chosen panel of the propeller

    #Loop for the number of blades
    for j in range (Var.Z_Blade_P):
        theta_blade = float(j*delta_theta)
        cos_theta = np.cos(theta_blade)
        sin_theta = np.sin(theta_blade)

        x_1 = x_10    # X value for the first point of the chosen panel of the chosen blade of the propeller
        y_1 = y_10*cos_theta - z_10*sin_theta
        # Y value for the first point of the chosen panel of the chosen blade of the propeller
        z_1 = z_10*cos_theta + y_10*sin_theta
        # Z value for the first point of the chosen panel of the chosen blade of the propeller

        # 4 sides of the panel
        for i in range(4):
            i_2 = i + 1 if i < 3 else 0

            x_2  = Grid_Points_P[N_Panel_P[n_pnl,i_2],0] # X value for the second point of the chosen panel
                                                         # of the chosen blade of the propeller
            y_20 = Grid_Points_P[N_Panel_P[n_pnl,i_2],1] # Y value for the second point of the chosen panel
                                                         # of the chosen blade of the propeller
            z_20 = Grid_Points_P[N_Panel_P[n_pnl,i_2],2] # Z value for the second point of the chosen panel
                                                         # of the chosen blade of the propeller

            y_2 = y_20 * cos_theta - z_20 * sin_theta
            # Y value for the second point of the chosen panel of the chosen blade of the propeller
            z_2 = z_20 * cos_theta + y_20 * sin_theta
            # Z value for the second point of the chosen panel of the chosen blade of the propeller

            if i == 1:
                x_1, y_1, z_1 = x_2, y_2, z_2
                continue
            if i == 3:
                x_1, y_1, z_1 = x_2, y_2, z_2
                continue
            else:
                U_x_0, U_y_0, U_z_0 = Biot_Savart_Propeller(1, x_1, y_1, z_1, x_2, y_2, z_2, px, py, pz)

                # Update induced velocities
                U_x = U_x + U_x_0
                U_y = U_y + U_y_0
                U_z = U_z + U_z_0

                x_1, y_1, z_1 = x_2, y_2, z_2

    return(U_x, U_y, U_z)
