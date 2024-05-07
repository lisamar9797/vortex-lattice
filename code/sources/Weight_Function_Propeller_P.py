"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine computes the weight function for the propeller,
involving the declaration of variables and arrays.
"""


import numpy as np
import sources.Variables as Var


def Weight_function_propeller():


    t_gp_P = np.loadtxt("output/Propeller_t_gp.txt", skiprows=1)
    # Rooftop parameter a
    a_roof = 0.8
    # 2 / Pi
    pi_inv = 2/np.pi
    # Domain limit of the distribution of circulation (rooftop)
    t_rest = 0.5 - a_roof
    # First Denominator
    t_slop = 2/(1 - a_roof*a_roof)
    # Second Denominator
    pcst = 2/(1 + a_roof)
    # Rooftop Coefficient
    cny1 = 1 - Var.cny

    GF_tot = [0.0]                         # Weight equation´s numerator
    GF_tmp = np.array([0.0]*(Var.Nch))     # Temporary Weight Function 2
    GW = np.array([0.0]*(Var.Nch))         # Temporary Weight Funtion 2
    G_Faux =np.array([0.0]*(Var.Nch))
    Weight_P = np.zeros((Var.Msp,Var.Nch))

    #Weight equation´s numerator
    for i in range(Var.Nch+1):
        i_1 = i-1
        t_gp_P_1 = 0.5 + t_gp_P[i]             # Numerator - Flat plate distribution (gamma)
        t_gp_P_2 = 0.5 - t_gp_P[i]             # Denominator - Flat plate distribution (gamma)
        Gamma_fp =pi_inv * np.sqrt(t_gp_P_1/t_gp_P_2) * Var.cny  # Flat plate distribution (gamma)
        Gamma_rt = pcst * cny1                                   # Rooftop distribution (gamma)

        if t_gp_P[i]<t_rest:                                     # Rooftop distribution (gamma)
            Gamma_rt = t_slop * t_gp_P_1 * cny1

        Gamma_Tot = Var.cny*Gamma_fp + cny1*Gamma_rt # This is the combination of flat plate distribution
                                                     # and rooftop distribution (gamma)

        GF_tmp[i_1] = Gamma_Tot*np.sqrt(t_gp_P_1*t_gp_P_2)
        GF_tot = GF_tot + GF_tmp[i_1]                # Weight equation´s denominator

    # Loop for the weight function (circulation of the panels)
    for i in range(Var.Nch-1,-1,-1):
        G_Faux += GF_tmp[i]/GF_tot
        GW[i] = G_Faux[i]
    #Loop used to write the weight function for all the panels# (it is the same for each spanwise level)
    for j in range (Var.Msp):
        for i in range (Var.Nch):
            Weight_P[j,i] = GW[i]   # Weight Function - j = spanwise level - i = chordwise level - Propeller

    with open("output/Propeller_Weight_Function.txt", "w") as file:
        file.write("Panel  Weight Function\n")
        for i in range(Var.Nch):
            file.write(f"{i:3d}    {Weight_P[0,i]:13.9f}\n")
    return(Weight_P)


Weight_P = Weight_function_propeller()
