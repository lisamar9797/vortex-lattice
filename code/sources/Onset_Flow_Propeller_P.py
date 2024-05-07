"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine computes the onset flow at the midpoints of the
sides of the propeller. The onset flow is assumed to be axi-symmetric and
independent of the longitudinal position, meaning it has only a radial
variation.
"""


import numpy as np
import sources.Variables as Var
from sources.Mid_Vect_Propeller_P import Mid_Vect_Propeller


def Onset_Flow_Propeller():
    I_P_Points_P = (Var.Msp*Var.Nch)
    r_R_P, X_P, Skew_P, Chord_P, Thick_P = np.loadtxt("input/grid.txt", unpack=True)
    U_0_P, U_R_P, U_T_P = np.loadtxt("input/onset.txt", unpack=True)
    U_0_P_Onset = np.zeros((Var.Msp * Var.Nch, 4, 3))   #Onset Flow (x) - Propeller
    U_T_P_Onset = np.zeros((Var.Msp * Var.Nch, 4, 3))   #Onset Flow (y) - Propeller
    U_R_P_Onset = np.zeros((Var.Msp * Var.Nch, 4, 3))   #Onset Flow (z) - Propeller
    V_Onset_P = np.zeros((Var.Msp * Var.Nch, 4, 3))

    for j in range(I_P_Points_P): # This loop selects the panel where the point px,py,pz is located
        for k in range (4):       # This loop selects, inside the panel, the sides where the point px,py,pz is located
            xx, xy, xz, xl, yl, zl = Mid_Vect_Propeller(j, k)

            #This subroutine is used to calculate the midpoint px,py,pz
            r_sid = np.sqrt(xy*xy + xz*xz)          # Radius of the points px,py,pz

            U_0_P_Onset= np.interp(r_sid,r_R_P,U_0_P)	 # Wake (Axial) in the midpoints (s)
            U_T_P_Onset= np.interp(r_sid,r_R_P,U_T_P)	 # Wake (Tangential) in the midpoints (s)
            U_R_P_Onset= np.interp(r_sid,r_R_P,U_R_P)    # Wake (Radial) in the grid midpoints (s)

            V_Onset_P[j,k,0] = - U_0_P_Onset		 # Onset Flow (x)
            V_Onset_P[j,k,1] = U_R_P_Onset*xy/r_sid - U_T_P_Onset*xz/r_sid + Var.Omega*xz  # Onset Flow (y)
            V_Onset_P[j,k,2] = U_R_P_Onset*xz/r_sid + U_T_P_Onset*xy/r_sid - Var.Omega*xy  # Onset Flow (z)

    with open("output/Propeller_Onset_Flow.txt", "w") as file:
        file.write("{:5s}{:>6s}{:12s}{:15s}\n".format("Point", "Ux", "Uy", "Uz"))
        file.write("{:<8s}{:<7s}\n".format("(Panel)", "(Side)"))

        for j in range(I_P_Points_P):
            for k in range(4):
                file.write("{:2d}{:4d}{:5s}{:13.9f}{:5s}{:13.9f}{:5s}{:13.9f}\n".format(
                    j, k, "", V_Onset_P[j, k, 0], "", V_Onset_P[j, k, 1], "", V_Onset_P[j, k, 2]))

    return V_Onset_P


V_Onset_P = Onset_Flow_Propeller()
