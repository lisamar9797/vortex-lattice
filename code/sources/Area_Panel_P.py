"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine calculates the area of a panel given its four
points.
"""


import numpy as np


def Area_Panel(x_1, y_1, z_1, x_2, y_2, z_2, x_3, y_3, z_3, x_4, y_4, z_4):


    s = 0                          # Initialization of the variable

    b_1 = x_4 - x_1                   # X value of the first vector of the panel (Point 1 and Point 4)
    b_2 = y_4 - y_1                   # Y value of the first vector of the panel (Point 1 and Point 4)
    b_3 = z_4 - z_1                   # Z value of the first vector of the panel (Point 1 and Point 4)

    e_1 = x_3 - x_1                   # First side (x)
    e_2 = y_3 - y_1                   # First side (y)
    e_3 = z_3 - z_1                   # First side (z)

    f_1 = x_2 - x_1                   # Second side (x)
    f_2 = y_2 - y_1                   # Second side (y)
    f_3 = z_2 - z_1                   # Second side (z)

    s_11 = f_2*b_3 - f_3*b_2          # X component of the first cross product
    s_12 = b_1*f_3 - f_1*b_3          # Y component of the first cross product
    s_13 = f_1*b_2 - f_2*b_1          # Z component of the first cross product

    s_21 = b_2*e_3 - b_3*e_2          # X component of the second cross product
    s_22 = e_1*b_3 - b_1*e_3          # Y component of the second cross product
    s_23 = b_1*e_2 - b_2*e_1          # Z component of the second cross product

    s = 0.5*(np.sqrt(s_11**2 + s_12**2 + s_13**2) + np.sqrt(s_21**2 + s_22**2 + s_23**2))  #Area of the panel

    return (s)
