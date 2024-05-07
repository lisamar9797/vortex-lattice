"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine calculates the advance ratio.
"""


import sources.Variables as Var
import numpy as np


def Advance_Ratio_J():

    r_R_P,X_P,Skew_P,Chord_P,Thick_P = np.loadtxt("input/grid.txt", unpack=True)
    U_0_P, U_R_P, U_T_P = np.loadtxt("input/onset.txt", unpack=True)

    n_0 = 50          # Number of intervals (It is used to find the "step length" for the composite Simpson's rule)
    n_1 = n_0         # Number of approximation values of the integral for the composite Simpson's rule
    h = (Var.Rad_P - Var.R_Hub_P)/n_0 # "step length"
    r_tmp = Var.R_Hub_P               # Initial value for the radius r

    Ua_tmp = 0          # Initialization of the approximation of the integral
    j = 1               # First value of j

    for i in range(n_1+1):                # Simpson's rule used in order to solve the integral
        j = - j                           #It used to have 2 or 4 in the composite Simpson's rule
        Simpson = 3 + float(j)            # This value is 2 or 4
        if(i == 0 or i == n_1):           # This if is used to have 1 as coefficient if we are considering the first
            Simpson = 1                   # or the last value of the integral

        Ux_tmp = np.interp(r_tmp,r_R_P,U_0_P)
        # Linear interpolation used to find the value of the axial velocity
        Ua_tmp = Ua_tmp + (Ux_tmp * r_tmp) * Simpson
        # Composite Simpson's rule
        r_tmp = r_tmp + h

    # Advance velocity
    U_adv = 2 * h * Ua_tmp /(3*(Var.Rad_P**2 - Var.R_Hub_P**2))
    # Advance ratio
    Advance_ratio =  (U_adv * np.pi) / (Var.Omega * Var.Rad_P)
    # Wake fraction
    w_eff = 1 - U_adv/Var.V_Ship

    if (w_eff < 1.0 - 10):
        w_eff = 0

    # Open the file for writing
    with open("output/Propeller_Hydrodynamic_Characteristics.txt", "w") as file:
        # Write the header line
        file.write("{:2s}{:16s}{:4s}{:13s}{:4s}{:13s}\n".format("", "Advance velocity", "", "Advance ratio", "", "Wake fraction"))
        file.write("{:3s}{:13.9f}{:5s}{:13.9f}{:4s}{:13.9f}\n".format("", U_adv, "", Advance_ratio, "", w_eff))
    return Advance_ratio


Advance_ratio = Advance_Ratio_J()
