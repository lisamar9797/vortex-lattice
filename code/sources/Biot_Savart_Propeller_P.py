"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine calculates the velocities Ux, Uy, Uz from a line
element (x1,y1,z1) to (x2,y2,z2) at the point (px,py,pz) using Biot-Savart's
law. The process applies to all blades of the propeller. Given the propeller's
symmetry, only the point on the reference blade is needed. The induced velocity
from the vortex is calculated with a unit circulation.
"""


import numpy as np


def Biot_Savart_Propeller(I_z, x_1, y_1, z_1, x_2, y_2, z_2, px, py, pz):

    U_x, U_y, U_z = 0.0, 0.0, 0.0

    d_theta = (2*np.pi)/float(I_z)  # Angle between the blades
    theta = 0.0                     # Angle for the first blade
    a_x = x_2 - x_1
    d_y = y_2 - y_1
    d_z = z_2 - z_1

    b_x = px - x_2
    c_x = px - x_1

    for i in range(I_z):
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)

        a_y = d_y * cos_theta - d_z * sin_theta
        a_z = d_z * cos_theta + d_y * sin_theta

        b_y = py - y_2 * cos_theta + z_2 * sin_theta
        b_z = pz - z_2 * cos_theta - y_2 * sin_theta

        c_y = py - y_1 * cos_theta + z_1 * sin_theta
        c_z = pz - z_1 * cos_theta - y_1 * sin_theta

        a_length = np.sqrt(a_x*a_x + a_y*a_y + a_z*a_z)  # Lenght a
        b_length = np.sqrt(b_x*b_x + b_y*b_y + b_z*b_z)  # Lenght b
        c_length = np.sqrt(c_x*c_x + c_y*c_y + c_z*c_z)  # Lenght c

        a_c = a_x*c_x + a_y*c_y + a_z*c_z  # Dot product a.c (e)
        a_b = a_x*b_x + a_y*b_y + a_z*b_z  # Dot product a.b (c-e)

        ac_x = a_y*c_z - a_z*c_y  # X component of the cross product a^c
        ac_y = a_z*c_x - a_x*c_z  # Y component of the cross product a^c
        ac_z = a_x*c_y - a_y*c_x  # Z component of the cross product a^c

        aclen2 = ac_x*ac_x + ac_y*ac_y + ac_z*ac_z
        aclen = np.sqrt(aclen2)  # Module of the cross product a^c

        # This if is used to check the distance between the selected point and
        # the side. If they are too close we have to skip it
        if a_length != 0 and (aclen / a_length) > 1*10**(-5):

            cstac = a_c/c_length  # e/c
            cstab = a_b/b_length  # a-e / b
            cstv = 1.0 / (4.0 * np.pi * aclen2)
            cstv1 = cstv*cstac - cstv*cstab

            U_x = U_x + ac_x * cstv1  # Induced Velocity (x)
            U_y = U_y + ac_y * cstv1  # Induced Velocity (y)
            U_z = U_z + ac_z * cstv1  # Induced Velocity (z)

        theta = theta + d_theta

    return (U_x, U_y, U_z)
