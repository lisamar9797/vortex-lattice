"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This function saves the old propeller pitch to evaluate the
residual for the pitch distribution.
"""


import numpy as np
import sources.Variables as Var


pitch_0 = np.zeros((Var.Msp+1, 1))


def pitch():
    pitch_0 = np.zeros((Var.Msp+1,1))
    Points_Trans_Wake_P = np.loadtxt("output/Propeller_Points_Trans_Wake.txt", skiprows= 1, usecols= (1,2,3))
    for i in range (Var.Msp+1):
        i_1 = i+i*(Var.N_P_L)
        pitch_0[i] = Points_Trans_Wake_P[i_1,2]
    return pitch_0
