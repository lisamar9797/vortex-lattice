"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine calculates the total velocities at the midpoints
of segments (Onset + Induced & Onset) for the reference blade of the propeller.
"""


from sources.Onset_Flow_Propeller_P import Onset_Flow_Propeller
from sources.Velocity_Total_No_Onset_Propeller_P import Velocity_Total_No_Onset_Propeller
import numpy as np
import sources.Variables as Var


def Velocity_Total_Propeller ():
    Gamma_TE_P = np.loadtxt("output/Propeller_Gamma_TE_P.txt")

    V_Ind_P, V_Tral_P = Velocity_Total_No_Onset_Propeller()
    V_Onset_P = Onset_Flow_Propeller()

    I_P_Points_P = (Var.Msp*Var.Nch)
    V_Tot_No_Onset_P = np.zeros((I_P_Points_P,4,3))
    V_Tot_P = np.zeros((I_P_Points_P,4,3))

    for i in range (I_P_Points_P):
        # This loop in used to select the panel where the point px,py,pz is located on the propeller

        # This loop in used to select the side of the panel where the point px,py,pz is located on the propeller
        for k in range(4):

            u_x = V_Onset_P[i,k,0]    # Temporary variable used to store the onset flow (x)
            u_y = V_Onset_P[i,k,1]    # Temporary variable used to store the onset flow (y)
            u_z = V_Onset_P[i,k,2]    # Temporary variable used to store the onset flow (z)

            u_xx = 0        # Initialization of the variable used to store the induced velocity
                            # of the propeller without the onset flow (x)
            u_yy = 0        # Initialization of the variable used to store the induced velocity
                            # of the propeller without the onset flow (y)
            u_zz = 0        # Initialization of the variable used to store the induced velocity
                            # of the propeller without the onset flow (z)

            # This loop is used to calculate the total velocity
            # at the midpoint of the panel of the propeller (Spanwise)
            for j in range (Var.Msp):
                u_x = u_x + Gamma_TE_P[j] * V_Ind_P[j,i,k,0]
                # Onset Flow + induced velocity in the panel i side k of the propeller (x)
                u_y = u_y + Gamma_TE_P[j] * V_Ind_P[j,i,k,1]
                # Onset Flow + induced velocity in the panel i side k of the propeller (y)
                u_z = u_z + Gamma_TE_P[j] * V_Ind_P[j,i,k,2]
                # Onset Flow + induced velocity in the panel i side k of the propeller (z)

                # Induced velocity in the panel i side k of the propeller
                u_xx = u_xx + Gamma_TE_P[j] * V_Ind_P[j,i,k,0]
                u_yy = u_yy + Gamma_TE_P[j] * V_Ind_P[j,i,k,1]
                u_zz = u_zz + Gamma_TE_P[j] * V_Ind_P[j,i,k,2]

            # Induced velocity of the propeller without the onset flow
            V_Tot_No_Onset_P[i,k,0]  = u_xx
            V_Tot_No_Onset_P[i,k,1]  = u_yy
            V_Tot_No_Onset_P[i,k,2]  = u_zz

            # Total induced velocity from the propeller in the panel i side k of the propeller
            V_Tot_P[i,k,0] = u_x
            V_Tot_P[i,k,1] = u_y
            V_Tot_P[i,k,2] = u_z


    # Open the file for Propeller_Velocity_Total_No_Onset
    with open("output/Propeller_Velocity_Total_No_Onset_V.txt", "w") as file:
        file.write("  Point                Ux            Uy              Uz\n")
        file.write("(Panel)  (Side)\n")

        for i in range(I_P_Points_P):
            for k in range(4):
                file.write(f"  {i:2d}    {k:4d}    {V_Tot_No_Onset_P[i, k, 0]:13.9f}   {V_Tot_No_Onset_P[i, k, 1]:13.9f}   {V_Tot_No_Onset_P[i, k, 2]:13.9f}\n")

    # Open the file for Propeller_Velocity_Trailing_Vortices
    with open("output/Propeller_Velocity_Total.txt", "w") as file:
        file.write("  Point                Ux            Uy              Uz\n")
        file.write("(Panel)  (Side)\n")

        for i in range(I_P_Points_P):
            for k in range(4):
                file.write(f"  {i:2d}    {k:4d}    {V_Tot_P[i,k,0]:13.9f}   {V_Tot_P[i, k, 1]:13.9f}   {V_Tot_P[i, k, 2]:13.9f}\n")

    return (V_Tot_P, V_Tot_No_Onset_P)


V_Tot_P, V_Tot_No_Onset_P = Velocity_Total_Propeller()
