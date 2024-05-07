"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine calculates the induced velocities (coefficient)
from a panel at the point (px,py,pz) for all blades of the propeller.

Parameters:
- n_pnl: Number of the panel inducing velocity at the point (px,py,pz)
- mpnl: Number of the panel containing the point (px,py,pz)
- msid: Number of the side containing the point (px,py,pz)
"""


from sources.Biot_Savart_Propeller_P import Biot_Savart_Propeller
import numpy as np
import sources.Variables as Var


def Panel_Induced_Velocity_Propeller(n_pnl, msid, mpnl, px, py, pz):
    Grid_Points_P = np.loadtxt("output/Propeller_Grid_Points.txt")
    N_Panel_P = np.loadtxt("output/Propeller_Numeration_Panel.txt", dtype='int')

    U_x, U_y, U_z = 0.0, 0.0, 0.0           # Initialization of the variable U_x, U_y, U_z

    delta_theta = 2*np.pi/float(Var.Z_Blade_P)

    x_10 = Grid_Points_P[N_Panel_P[n_pnl,0],0]  # X value for the first point of the chosen panel of the propeller
    y_10 = Grid_Points_P[N_Panel_P[n_pnl,0],1]  # Y value for the first point of the chosen panel of the propeller
    z_10 = Grid_Points_P[N_Panel_P[n_pnl,0],2]  # Z value for the first point of the chosen panel of the propeller

    for j in range(Var.Z_Blade_P):    #Loop for the number of blades
        theta_blade = (j) *delta_theta

        cos_theta = np.cos(theta_blade)
        sin_theta = np.sin(theta_blade)

        x_1 = x_10
        y_1 = y_10*cos_theta - z_10*sin_theta
        z_1 = z_10*cos_theta + y_10*sin_theta
        # X value for the first point of the chosen panel of the chosen blade of the propeller
        # Y value for the first point of the chosen panel of the chosen blade of the propeller
        # Z value for the first point of the chosen panel of the chosen blade of the propeller

        # 4 sides of the panel
        for i in range(4):

            i_2 = i + 1 if i < 3 else 0
            x_2  = Grid_Points_P[N_Panel_P[n_pnl,i_2],0]
            # X value for the second point of the chosen panel of the chosen blade of the propeller
            y_20 = Grid_Points_P[N_Panel_P[n_pnl,i_2],1]
            # Y value for the second point of the chosen panel of the chosen blade of the propeller
            z_20 = Grid_Points_P[N_Panel_P[n_pnl,i_2],2]
            # Z value for the second point of the chosen panel of the chosen blade of the propeller

            y_2 = y_20 * cos_theta - z_20 * sin_theta
            # Y value for the second point of the chosen panel of the chosen blade of the propeller
            z_2 = z_20 * cos_theta + y_20 * sin_theta
            # Z value for the second point of the chosen panel of the chosen blade of the propeller

            if n_pnl == mpnl and (j == 0) and i == msid:
                x_1 = x_2               # The second point becomes the first point (x)
                y_1 = y_2               # The second point becomes the first point (y)
                z_1 = z_2               # The second point becomes the first point (z)
            else:
                U_x_0, U_y_0, U_z_0 = Biot_Savart_Propeller(1, x_1, y_1, z_1, x_2, y_2, z_2, px, py, pz)
                # Induced velocity of that side of the panel of the propeller
                # (I_z = 1 because I have al ready created a loop)
                # Update induced velocity

                U_x = U_x + U_x_0
                U_y = U_y + U_y_0
                U_z = U_z + U_z_0

                x_1 = x_2               # The second point becomes the first point (x)
                y_1 = y_2               # The second point becomes the first point (y)
                z_1 = z_2               # The second point becomes the first point (z)

    return(U_x, U_y, U_z)
