"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine calculates the midpoint coordinates
(mid_x,mid,mid_z) and the vector for the panel side
(vector_x,vector_y,vector_z) of the reference blade of the propeller.

Parameters
 - n_pnl: number of the panel
 - n_side: number of the side
"""


import numpy as np
import sources.Variables as Var


def Mid_Vect_Propeller(n_pnl,n_side):

    Grid_Points_P = np.loadtxt("output/Propeller_Grid_Points.txt")
    N_Panel_P = np.loadtxt("output/Propeller_Numeration_Panel.txt",dtype='int')
    j1 = n_side
    j2 = (n_side+1)

    # Special case for n_side = 3
    if n_side == 3:
        j2 = 0

    mid_x = 0.5* (Grid_Points_P[N_Panel_P[n_pnl,j2],0] + Grid_Points_P[N_Panel_P[n_pnl,j1],0])		# Midpoint (x)
    mid_y = 0.5* (Grid_Points_P[N_Panel_P[n_pnl,j2],1] + Grid_Points_P[N_Panel_P[n_pnl,j1],1])		# Midpoint (y)
    mid_z = 0.5* (Grid_Points_P[N_Panel_P[n_pnl,j2],2] + Grid_Points_P[N_Panel_P[n_pnl,j1],2])      # Midpoint (z)

    vector_x = Grid_Points_P[N_Panel_P[n_pnl,j2],0] - Grid_Points_P[N_Panel_P[n_pnl,j1],0]		#Vector (x)
    vector_y = Grid_Points_P[N_Panel_P[n_pnl,j2],1] - Grid_Points_P[N_Panel_P[n_pnl,j1],1]		#Vector (y)
    vector_z = Grid_Points_P[N_Panel_P[n_pnl,j2],2] - Grid_Points_P[N_Panel_P[n_pnl,j1],2]		#Vector (z)

    return mid_x,mid_y,mid_z,vector_x,vector_y,vector_z
