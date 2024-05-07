"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine is tasked with creating the helix geometry.
"""


def Helix(x1, x2, y1, y2, dx):
    delx = x2 - x1
    a = (y2-y1-delx*dx)/delx/delx
    b = dx-2*a*x1
    c = y1 + a*x1*x1 - dx*x1
    return (a, b, c)
