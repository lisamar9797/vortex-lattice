"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine calculates the Ci-function used in 'De_Jong' by
rational approximations.
"""


import numpy as np


def Ci(xbar):


    f = (xbar**8+38.027264*xbar**6+265.187033*xbar**4+335.677320*xbar**2+38.102495)/(
    xbar**8+40.021433*xbar**6+322.624911*xbar**4+570.236280*xbar**2+157.105423)/xbar

    g = (xbar**8+42.242855*xbar**6+302.757865*xbar**4+352.018498*xbar**2+21.821899)/(
    xbar**8+48.196927*xbar**6+482.485984*xbar**4+1114.978885*xbar**2+449.690326)/xbar**2

    Cires=f*np.sin(xbar)-g*np.cos(xbar)

    return(Cires)
