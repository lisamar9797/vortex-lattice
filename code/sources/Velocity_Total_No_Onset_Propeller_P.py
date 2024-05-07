"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine computes the induced velocities at the midpoints
of segments (coefficient) emanating from the horseshoe vortex across all blades
of the propeller. Furthermore, it calculates the total "velocity matrix" at
these midpoints, excluding the onset flow.
"""


import numpy as np
import sources.Variables as Var
from sources.Induced_Grid_Propeller_P import Induced_Grid_Propeller
from sources.Mid_Vect_Propeller_P import Mid_Vect_Propeller
from sources.Biot_Savart_Propeller_P import Biot_Savart_Propeller
from sources.Trailing_Vortices_Propeller_P import Trailing_Vortices_Propeller


def Velocity_Total_No_Onset_Propeller():

    Horseshoe_P = np.loadtxt("output/Propeller_Horseshoe.txt", dtype='int')
    Grid_Points_P = np.loadtxt("output/Propeller_Grid_Points.txt")
    N_Bound_Vortex_P = np.loadtxt("output/Propeller_N_Bound_Vortex.txt", dtype='int')
    N_Bound_Vortex_P = N_Bound_Vortex_P.reshape((Var.Msp+1, 1))

    I_P_Points_P = (Var.Msp*Var.Nch)
    V_Grid_P = Induced_Grid_Propeller()
    V_Tral_P = np.zeros((Var.Msp+1, I_P_Points_P, 4,3))
    V_Ind_P = np.zeros((Var.Msp, I_P_Points_P, 4,3))

    # HORSESHOE VORTEX

    for i in range (I_P_Points_P):
    # I used this loop in order to select the panel where the point xx,xy,xz is located

        N_P_V = Horseshoe_P[0,0]
        # First point of the first trailing vortex selected (Segment)

        for k in range (4):
        # I used this loop in order to select the side where the point xx,xy,xz
        #is located and to calculate the induced velocity
    	# from the transition wake and from the semi-infinite helicoidal vortex for the selected trailing vortex (N_P_V)
            xx,xy,xz,v_xx,v_xy,v_xz = Mid_Vect_Propeller(i,k)   # This subroutine is used to calculate the midpoint xx,xy,xz

            U_x1,U_y1,U_z1 = Trailing_Vortices_Propeller (N_P_V,xx,xy,xz)
            # Induced velocity from the transition wake and from the semi-infinite helicoidal vortex (First)

            V_Tral_P[0,i,k,0] = U_x1
            V_Tral_P[0,i,k,1] = U_y1
            V_Tral_P[0,i,k,2] = U_z1

        x_1 = Grid_Points_P[N_Bound_Vortex_P[N_P_V,0],0] # X coordinate of the second point of the segment of the first trailing vortex
        y_1 = Grid_Points_P[N_Bound_Vortex_P[N_P_V,0],1] # Y coordinate of the second point of the segment of the first trailing vortex
        z_1 = Grid_Points_P[N_Bound_Vortex_P[N_P_V,0],2] # Z coordinate of the second point of the segment of the first trailing vortex

        for j in range(Var.Msp):						# This loop is used to select the horseshoe vortex

            j_2 = j+1
            N_P_V_2 = Horseshoe_P[j,1]    # Second point of the trailing vortex selected (Segment)

            x_2 = Grid_Points_P[N_Bound_Vortex_P[N_P_V_2,0],0]	# X coordinate of the second point of the segment of the trailing vortex
            y_2 = Grid_Points_P[N_Bound_Vortex_P[N_P_V_2,0],1]	# Y coordinate of the second point of the segment of the trailing vortex
            z_2 = Grid_Points_P[N_Bound_Vortex_P[N_P_V_2,0],2]	# Z coordinate of the second point of the segment of the trailing vortex

            for k in range (4):
                xx,xy,xz,v_xx,v_xy,v_xz= Mid_Vect_Propeller(i,k)

                U_x2,U_y2,U_z2 = Trailing_Vortices_Propeller(N_P_V_2,xx,xy,xz)
                # Induced velocity from the transition wake and from the semi-infinite helicoidal vortex (Second)

                U_xs,U_ys,U_zs = Biot_Savart_Propeller(Var.Z_Blade_P,x_1,y_1,z_1,x_2,y_2,z_2,xx,xy,xz)
                # Induced velocity from the bound vortex selected

                V_Tral_P[j,i,k,0] = V_Tral_P[j,i,k,0] - U_x2 + U_xs  # X velocity induced from the horseshoe vortex
                V_Tral_P[j,i,k,1] = V_Tral_P[j,i,k,1] - U_y2 + U_ys  # Y velocity induced from the horseshoe vortex
                V_Tral_P[j,i,k,2] = V_Tral_P[j,i,k,2] - U_z2 + U_zs	 # Z velocity induced from the horseshoe vortex


                V_Tral_P[j_2,i,k,0] = U_x2     # I need this value for the next i loop (U_x2 will be U_x1 for the next horseshoe vortex)
                V_Tral_P[j_2,i,k,1] = U_y2	   # I need this value for the next i loop (U_y2 will be U_z1 for the next horseshoe vortex)
                V_Tral_P[j_2,i,k,2] = U_z2	   # I need this value for the next i loop(U_y2 will be U_z1 for the next horseshoe vortex)

            x_1 = x_2           # This is used in order to have the first point of the next bound vortex (x)
            y_1 = y_2           # This is used in order to have the first point of the next bound vortex (y)
            z_1 = z_2           # This is used in order to have the first point of the next bound vortex (z)

    # VELOCITY MATRIX

    for j in range(Var.Msp):
        for i in range(I_P_Points_P):
            for k in range (4):
                V_Ind_P[j,i,k,0] = V_Grid_P [j,i,k,0] + V_Tral_P[j,i,k,0]  # Total induced velocity without the onset flow (x)
                V_Ind_P[j,i,k,1] = V_Grid_P [j,i,k,1] + V_Tral_P[j,i,k,1]  # Total induced velocity without the onset flow (y)
                V_Ind_P[j,i,k,2] = V_Grid_P [j,i,k,2] + V_Tral_P[j,i,k,2]# Total induced velocity without the onset flow (y)

    # Open the file for Propeller_Velocity_Total_No_Onset
    with open("output/Propeller_Velocity_Total_No_Onset.txt", "w") as file:
        file.write("    Point         Spanwise       Ux              Uy               Uz\n")
        file.write("  (Panel)  (Side)\n")

        for i in range(I_P_Points_P):
            for k in range(4):
                for j in range(Var.Msp):
                    file.write(f"    {i:2d}    {k:4d}     {j:4d}    {V_Ind_P[j, i, k, 0]:13.9f}    {V_Ind_P[j, i, k, 1]:13.9f}   {V_Ind_P[j, i, k, 2]:13.9f}\n")

    # Open the file for Propeller_Velocity_Trailing_Vortices
    with open("output/Propeller_Velocity_Trailing_Vortices.txt", "w") as file:
        file.write("    Point         Spanwise       Ux              Uy               Uz\n")
        file.write("(Panel)   (Side)\n")

        for i in range(I_P_Points_P):
            for k in range(4):
                for j in range(Var.Msp):
                    file.write(f"    {i:2d}    {k:4d}     {j:4d}    {V_Tral_P[j, i, k,0]:13.9f}    {V_Tral_P[j, i, k, 1]:13.9f}   {V_Tral_P[j, i, k, 2]:13.9f}\n")

    return V_Ind_P, V_Tral_P


V_Ind_P, V_Tral_P = Velocity_Total_No_Onset_Propeller()
