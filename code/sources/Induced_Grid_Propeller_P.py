"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine computes the induced velocities in the midpoints
of the segments (coefficient) from the entire grid of the propeller.
"""


import numpy as np
import sources.Variables as Var
from sources.Weight_Function_Propeller_P import Weight_function_propeller
from sources.Mid_Vect_Propeller_P import Mid_Vect_Propeller
from sources.Panel_Induced_Velocity_Propeller_P import Panel_Induced_Velocity_Propeller


def Induced_Grid_Propeller():
    Weight_P = Weight_function_propeller()
    I_P_Points_P = (Var.Msp*Var.Nch)

    V_Grid_P = np.zeros((Var.Msp,I_P_Points_P, 4,3))  # Induced velocity from the entire grid
    n_plaux = np.array(([0.0]*(Var.Msp)), dtype = int)
    npl = np.array(([0.0]*(Var.Msp)),dtype = int)

    for i in range (I_P_Points_P): # This loop selects the panel where the point px,py,pz is located

        for k in range (4): # This loop selects, inside the panel, the side where the point px,py,pz is located
            px,py,pz,v_px,v_py,v_pz = Mid_Vect_Propeller(i, k) # This subroutine is used to calculate the midpoint px,py,pz

            for j in range(Var.Msp):
                U_x = 0                 # Initialization of the variable U_x
                U_y = 0                 # Initialization of the variable U_y
                U_z = 0                 # Initialization of the variable U_z
                n_plaux = (Var.Nch)*(j)

                for h in range(Var.Nch): # This loop selects the panel (Chordwise) that induces velocity

                    npl = h + n_plaux
                    qx_pnl, qy_pnl, qz_pnl = Panel_Induced_Velocity_Propeller(npl,k,i,px,py,pz)

                    U_x += Weight_P[j,h] * qx_pnl   # Temporary induced velocity in the point px,py,pz
                    U_y += Weight_P[j,h] * qy_pnl   # due to the chordwise ring j (x),(y),(z)
                    U_z += Weight_P[j,h] * qz_pnl

                    V_Grid_P [j,i,k,0] = U_x         # Induced velocity in the point px,py,pz
                    V_Grid_P [j,i,k,1] = U_y         # due to the chordwise ring j (x),(y),(z)
                    V_Grid_P [j,i,k,2] = U_z

    with open("output/Propeller_Velocity_Grid.txt", "w") as file:
        file.write("{:>5s} {:>8s} {:>2s} {:>15s} {:>15s} {:>2s}\n".format("Point","Spanwise","Ux","Uy","Uz", ""))
        file.write("{:>8s} {:>7s}\n".format("(Panel)", "(Side)"))

        for i in range(I_P_Points_P):
            for k in range(4):
                for j in range(Var.Msp):
                    data_format = "{:>2d} {:>4d} {:>4d} {:>13.9f} {:>13.9f} {:>13.9f}\n"
                    file.write(data_format.format(i, k, j, V_Grid_P[j, i, k, 0], V_Grid_P[j, i, k, 1], V_Grid_P[j, i , k, 2]))
    return V_Grid_P


V_Grid_P = Induced_Grid_Propeller()
