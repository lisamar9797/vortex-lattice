"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine calculates the normal vector for a panel
"""


import numpy as np


def Normal_Vector(x_1, y_1, z_1, x_2, y_2, z_2, x_3, y_3, z_3, x_4, y_4, z_4):
    a_1 = x_2 - x_3
    # X value of the first vector of the panel (Point 2 and Point 3)
    a_2 = y_2 - y_3
    # Y value of the first vector of the panel (Point 2 and Point 3)
    a_3 = z_2 - z_3
    # Z value of the first vector of the panel (Point 2 and Point 3)

    b_1 = x_4 - x_1
    # X value of the first vector of the panel (Point 1 and Point 4)
    b_2 = y_4 - y_1
    # Y value of the first vector of the panel (Point 1 and Point 4)
    b_3 = z_4 - z_1
    # Z value of the first vector of the panel (Point 1 and Point 4)

    x = a_2 * b_3 - a_3 * b_2  # X component of the cross product
    y = b_1 * a_3 - a_1 * b_3  # Y component of the cross product
    z = a_1 * b_2 - a_2 * b_1  # Z component of the cross product

    Norm = np.sqrt(x**2 + y**2 + z**2)  # Norm of the vector

    vector_x = x / Norm  # X component of the normal vector
    vector_y = y / Norm  # Y component of the normal vector
    vector_z = z / Norm  # Z component of the normal vector

    return (vector_x, vector_y, vector_z)
