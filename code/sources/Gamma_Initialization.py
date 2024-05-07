"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This section is dedicated to the initialization of the variable
Gamma_TE_P.
"""


import numpy as np
import sources.Variables as Var


def Gamma_It():

    Gamma_TE_P = np.zeros((Var.Msp+1))

    for j in range(Var.Msp+1):
        Gamma_TE_P[j] = 0.0
    Gamma_TE_P[Var.Msp] = -1      # lambda (t-1) initial
    with open("output/Propeller_Gamma_TE_P.txt", "w") as file:
        for i in range(Var.Msp+1):
            file.write(f"{Gamma_TE_P[i]:13.9f}\n")
    return Gamma_TE_P


Gamma_TE_P = Gamma_It()
