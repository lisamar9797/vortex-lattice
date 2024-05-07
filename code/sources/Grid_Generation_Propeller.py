"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine is responsible for creating the initial grid for
the reference blade of the propeller, including Control Points & Grid Points.
Given the propeller's symmetry, generating the grid for the reference blade
alone suffices. Additionally, the subroutine undertakes the numbering of panels
and horseshoe vortices. Notably, the grid is aligned with the onset flow during
this phase.
"""


import math
import numpy as np
import sources.Variables as Var
import pandas as pd


def Grid_Generation_Propeller():

    r_R_P, X_P, Skew_P, Chord_P, Thick_P = np.loadtxt("input/grid.txt", unpack=True)
    U_0_P, U_R_P, U_T_P = np.loadtxt("input/onset.txt", unpack=True)

    """ MID-CHORD LINE """

    Midchord_line_P = np.zeros((Var.N_Iter,3))  # Create Matrix 3xN_Iter
    # Initial point (x, y, z)
    Midchord_line_P[0,0] = 0.0       # s line (x) - It begins at the hub-center but the first point is at the top of the hub
    Midchord_line_P[0,1] = 0.0       # s line (y) - It begins at the hub-center but the first point is at the top of the hub
    Midchord_line_P[0,2] = r_R_P[0]  # s line (z) assuming that r_R_P contain the initial z values

    # Cartesian coordinates for the blade surface
    for i in range(1, Var.N_Iter):
        Midchord_line_P[i, 0] = X_P[i]                            # x
        Midchord_line_P[i, 1] = -r_R_P[i] * math.sin(Skew_P[i])   # y
        Midchord_line_P[i, 2] = r_R_P[i] * math.cos(Skew_P[i])    # z

    s_tip = 0.0
    S_Distr_P = [0.0] * Var.N_Iter
    S_Distr_P[0] = math.sqrt(Midchord_line_P[0,1]**2 + Midchord_line_P[0,2]**2) # First value of the midchord line
    s_Hub_P = S_Distr_P[0]

    for i in range( Var.N_Iter - 1):  # This loop is used to find the length of the s line
        b = abs(Midchord_line_P[i+1,1])-abs(Midchord_line_P[i,1])       # Y Distance
        c = abs(Midchord_line_P[i+1,2])-abs(Midchord_line_P[i,2])       #Z Distance
        Prov = math.sqrt(b**2+c**2)
        S_Distr_P[i+1] = S_Distr_P[i] + Prov  # This is used to find the distribution of s, which is always costant

    s_tip = S_Distr_P[Var.N_Iter-1]

    data = np.column_stack([S_Distr_P, r_R_P])
    np.savetxt("output/Propeller_S_Distr.txt", data, fmt=['%13.9f','%13.9f'], delimiter= ' ', header = '    S_Distr        Radius')

    """ t FUNCTION     """

    t_gp_P=np.array([0.0]*(Var.Nch+1))    #Initialise the variable

    for i in range(Var.Nch + 1):
        t_gp_P[0] = -0.5
        t_gp_P[i] = -0.5 * np.cos(float(i+1-1.5)*3.14159274/float(Var.Nch))
        # (t) Grid points (always the same - it depends on Nch) - Cosine

    t_cp_P = np.array([0.0]*(Var.Nch),dtype=np.float64)  #Initialise the variable
    for i in range (Var.Nch):
        t_cp_P[i]= 0.5*(t_gp_P[i+1]+t_gp_P[i])  # (t) Control points (always the same - it depends on Nch) - Cosine


    data_t = np.column_stack([t_gp_P])
    np.savetxt("output/Propeller_t_gp.txt", data_t, fmt=['%13.9f'], delimiter= ' ', header = 't_gp')
    data_t = np.column_stack([t_cp_P])
    np.savetxt("output/Propeller_t_cp.txt", data_t, fmt=['%13.9f'], delimiter= ' ', header = 't_cp')

    """  s FUNCTION     """

    s_gp_P = np.array([0.0]*(Var.Msp+1),dtype=np.float64)       # (s) Grid points (always the same - it depends on Msp)
    for i in range(Var.Msp+1):
        aa = (i+1)*4.0 - 3.0 #Start with point after 0
        bb = 4.0*float(Var.Msp) + 2.0
        s_gp_P[i] = ((aa/bb)*(s_tip - s_Hub_P))+s_Hub_P

    s_cp_P = np.array([0.0]*(Var.Msp),dtype=np.float64)         # (s) Control points (always the same - it depends on Msp)
    for i in range(0,Var.Msp):
        s_cp_P[i] = 0.5 * (s_gp_P[i] + s_gp_P[i+1])

    data_s = np.column_stack([s_gp_P])
    np.savetxt("output/Propeller_s_gp.txt", data_s, fmt=['%13.9f'], delimiter= ' ', header = 's_gp')

    data_s = np.column_stack([s_cp_P])
    np.savetxt("output/Propeller_s_cp.txt", data_s, fmt=['%13.9f'], delimiter= ' ', header = 's_cp')

    """  GRID POINTS MATRIX - CALCULATION OF BETA(S),CHORD(S),SKEW(S) AND RAKE(S)  """

    #Initialise the variables
    Radius_gp_P = np.zeros(Var.Msp + 1)
    Chord_P_gp = np.zeros(Var.Msp + 1)
    Rake_P_gp = np.zeros(Var.Msp + 1)
    Skew_P_gp = np.zeros(Var.Msp + 1)
    sin_b = np.zeros(Var.Msp + 1)
    cos_b = np.zeros(Var.Msp + 1)
    Grid_Points_P = np.zeros(((Var.Msp + 1) * (Var.Nch + 1), 3))
    Theta_gp_P = np.zeros(Var.Nch + 1)

    # S Loop
    for i in range( Var.Msp +1):
        ipl = ((i+1) * (Var.Nch + 1)) - (Var.Nch + 1)

        Radius_gp_P[i] = np.interp(s_gp_P[i],S_Distr_P, r_R_P)  # Value of the radius in the grid points (s)
        U_0_P_gp = np.interp(s_gp_P[i],S_Distr_P, U_0_P)        # Wake (Axial) in the grid points (s)
        U_T_P_gp = np.interp(s_gp_P[i],S_Distr_P, U_T_P)        # Wake (Tangential) in the grid points (s)
        Chord_P_gp[i] = np.interp(s_gp_P[i],S_Distr_P, Chord_P) # Value of the chord in the grid points (s)
        Rake_P_gp[i] = np.interp(s_gp_P[i],S_Distr_P, X_P)      # Value of the rake in the grid points (s)
        Skew_P_gp[i] = np.interp(s_gp_P[i],S_Distr_P, Skew_P)   # Value of the skew in the grid points (s)

        V_tang = Var.Omega * Radius_gp_P[i] - U_T_P_gp  # Tangential velocity
        V_rel = np.sqrt(V_tang**2 + U_0_P_gp**2)        # Tangential velocity
        sin_b[i] = U_0_P_gp / V_rel                     # Sine (beta)
        cos_b[i] = V_tang / V_rel                       # Cosine (beta)

        # t Loop
        for j in range(Var.Nch+1):
            npl = (j) + ipl # Second counter used to order the Grid Points Matrix
            Theta_gp_P[j] = -Skew_P_gp[i] + (t_gp_P[j] * Chord_P_gp[i] * cos_b[i]) / Radius_gp_P[i]

            Grid_Points_P[npl, 0] = Rake_P_gp[i] + Chord_P_gp[i] * sin_b[i] * t_gp_P[j]   # X(s,t)
            Grid_Points_P[npl, 1] = -Radius_gp_P[i] * math.sin(Theta_gp_P[j])             # Y(s,t)
            Grid_Points_P[npl, 2] = Radius_gp_P[i] * math.cos(Theta_gp_P[j])              # Z(s,t)

    with open('output/Propeller_Grid_Points_Old.txt', 'w') as file:
        for i in range((Var.Nch + 1) * (Var.Msp + 1)):
            file.write(f"{Grid_Points_P[i, 0]:.9f}      {Grid_Points_P[i, 1]:.9f}      {Grid_Points_P[i, 2]:.9f}\n")
    with open('output/Propeller_Grid_Points.txt', 'w') as file:
        for i in range((Var.Nch + 1) * (Var.Msp + 1)):
            file.write(f"{Grid_Points_P[i, 0]:.9f}      {Grid_Points_P[i, 1]:.9f}      {Grid_Points_P[i, 2]:.9f}\n")
    with open('output/Propeller_Grid_Points_geom.txt', 'w') as file:
        for i in range((Var.Msp + 1)):
            file.write(f"{Radius_gp_P[i]:.9f}      {Chord_P_gp[i]:.9f}      {Rake_P_gp[i]:.9f}    {Skew_P_gp[i]:.9f}\n")

    """  CONTROL POINTS MATRIX - CALCULATION OF BETA(S(R)),CHORD(S),SKEW(S) AND RAKE(S)  """

    #Initialise the variable
    Radius_cp_P = np.zeros(Var.Msp + 1)
    Chord_P_cp = np.zeros(Var.Msp + 1)
    Rake_P_cp = np.zeros(Var.Msp + 1)
    Skew_P_cp = np.zeros(Var.Msp + 1)
    Theta_cp_P = np.zeros(Var.Nch + 1)
    Control_Points_P = np.zeros(((Var.Msp) * (Var.Nch), 3))

    for i in range(Var.Msp):
        Radius_cp_P[i] = 0.5*(Radius_gp_P[i]+Radius_gp_P[i+1]) # Value of the radius in the control point (s)
        U_0_P_cp = np.interp(s_cp_P[i],S_Distr_P, U_0_P)       # Wake (Axial) in the control points (s)
        U_T_P_cp = np.interp(s_cp_P[i],S_Distr_P, U_T_P)       #Wake (Tangential) in the control points (s)
        ipl = [0.0]
        ipl = ((i+1)*(Var.Nch))-(Var.Nch)              # Counter used to order the Control Points Matrix

        Chord_P_cp [i] = np.interp(s_cp_P[i],S_Distr_P, Chord_P)
        Rake_P_cp [i] = np.interp(s_cp_P[i],S_Distr_P, X_P)
        Skew_P_cp [i] = np.interp(s_cp_P[i],S_Distr_P, Skew_P)
        # Value of the chord in the control point (s)
        # Value of the rake in the control point (s)
        # Value of the skew in the control point (s)


        V_tang = Var.Omega * Radius_cp_P[i] - U_T_P_cp        # Tangential velocity
        V_rel = math.sqrt(V_tang**2 + U_0_P_cp**2)            # Relative velocity
        sin_b[i] = U_0_P_cp/V_rel                             # Sine (beta)
        cos_b[i] = V_tang/V_rel							      # Cosine (beta)

        #t loop
        for j in range(Var.Nch):
            npl = j+(ipl) # Second counter used to order the Control Points Matrix
            Theta_cp_P[j] = -Skew_P_cp[i] + (t_cp_P[j] * Chord_P_cp[i] * cos_b[i]) / Radius_cp_P[i]
            Control_Points_P[npl, 0] = Rake_P_cp[i] + Chord_P_cp[i] * sin_b[i] * t_cp_P[j]    # X(s,t)
            Control_Points_P[npl, 1] = -Radius_cp_P[i] * math.sin(Theta_cp_P[j])              # Y(s,t)
            Control_Points_P[npl, 2] = Radius_cp_P[i] * math.cos(Theta_cp_P[j])               # Z(s,t)

    with open('output/Propeller_Control_Points_Old.txt', 'w') as file:
        for i in range((Var.Nch) * (Var.Msp)):
            file.write(f"{Control_Points_P[i, 0]:.9f}      {Control_Points_P[i, 1]:.9f}      {Control_Points_P[i, 2]:.9f}\n")

    with open('output/Propeller_Control_Points.txt', 'w') as file:
        for i in range((Var.Nch) * (Var.Msp)):
            file.write(f"{Control_Points_P[i, 0]:.9f}      {Control_Points_P[i, 1]:.9f}      {Control_Points_P[i, 2]:.9f}\n")

    with open('output/Propeller_Control_Points_geom.txt', 'w') as file:
        for i in range((Var.Msp+1)):
            file.write(f"{Chord_P_cp [i]:.9f}      {Rake_P_cp [i]:.9f}      {Skew_P_cp [i]:.9f}\n")

    """  NUMERATION OF THE PANEL AND THE SIDE  """

    N_Panel_P = np.array([[0,1, Var.Nch + 2, Var.Nch +1]])     # Initialize the first panel

    t = 0
    for j in range(Var.Msp):
        for i in range(1,Var.Nch):

            t += 1
            t2 = t - 1

            N_Panel = N_Panel_P[t2] + 1
            N_Panel_P = np.append(N_Panel_P, [N_Panel], axis=0)

        if j != Var.Msp-1:
            t += 1
            t1 = t - 1
            N_Panel = N_Panel_P[t1] + 2
            N_Panel_P = np.append(N_Panel_P, [N_Panel], axis=0)

    with open ("output/Propeller_Numeration_Panel.txt","w") as file:
        for i in range((Var.Nch*Var.Msp)):
            file.write(f"{N_Panel_P[i, 0]:2d}   {N_Panel_P[i, 1]:2d}   {N_Panel_P[i, 2]:2d}   {N_Panel_P[i, 3]:2d}\n")

    """  COORDINATES FOR THE BOUND VORTICES (T.E. SIDE)  """

    # It is a matrix with the grid points at the T.E.
    N_Bound_Vortex_P = np.zeros((Var.Msp + 1, 1), dtype=int)
    for i in range (Var.Msp+1):
        N_Bound_Vortex_P[i] = i+i*(Var.Nch)

    np.savetxt("output/Propeller_N_Bound_Vortex.txt", N_Bound_Vortex_P, fmt="%3d")

    """  HORSESHOE VORTEX MATRIX     """

    Horseshoe_P = np.zeros(((Var.Msp),4),dtype=int)
    for i in range(Var.Msp):
        Horseshoe_P[i,0] = i                       # This value is used to access to the N_Bound_Vortex_P
                                                   # matrix in order to select the Bound vortex (0-6)(6-12)..
        Horseshoe_P[i,1] = i+1 					   # This value is used to access to the N_Bound_Vortex_P
                                                   # matrix in order to select the Bound vortex (0-6)(6-12)..
        Horseshoe_P[i,2] =((i-1)+1)*Var.Nch	       # Number of the T.E panel (0-5-10..)
        Horseshoe_P[i,3] = i                       # Number of the horseshoe vortex (0-1-2-3..)

    with open("output/Propeller_Horseshoe.txt","w") as file:
        for i in range(Var.Msp):
            file.write(f"   {Horseshoe_P[i,0]:2d}   {Horseshoe_P[i,1]:2d}   {Horseshoe_P[i,2]:2d}    {Horseshoe_P[i,3]}\n")

    """  COORDINATES FOR THE TRANSITION WAKE (STRAIGHT LINE VORTICES)  """

    # Initialize the variable
    Points_Trans_Wake_P = np.zeros((((Var.N_P_L+1)*(Var.Msp+1)),3))

    for i in range(Var.Msp+1):
        i_1 = i+i*(Var.N_P_L)
        x_trans_wake = Grid_Points_P[N_Bound_Vortex_P[i,0],0]
        y_trans_wake = Grid_Points_P[N_Bound_Vortex_P[i,0],1]
        z_trans_wake = Grid_Points_P[N_Bound_Vortex_P[i,0],2]
        # X value for the first point of the transition wake - T.E.
        # Y value for the first point of the transition wake - T.E.
        # Z value for the first point of the transition wake - T.E.

        r_trans_wake = np.sqrt(y_trans_wake**2 + z_trans_wake**2) # Radius at the T.E.
        U_0_P_trans_wake = np.interp(r_trans_wake, r_R_P, U_0_P)  # Wake (Axial) in the transition wake (s)
        U_T_P_trans_wake = np.interp(r_trans_wake, r_R_P, U_T_P)  # Wake (Tangential) in the transition wake (s)


        V_tang = Var.Omega*r_trans_wake - U_T_P_trans_wake      # Tangential velocity
        pitch_trans_wake = (2*np.pi*r_trans_wake*U_0_P_trans_wake)/V_tang  # Pitch at the T.E. (It has only a radial variation)
        Points_Trans_Wake_P[i_1,0] = x_trans_wake         # Grid points for the transition wake (x) - T.E
        Points_Trans_Wake_P[i_1,1] = r_trans_wake         # Grid points for the transition wake (radius) - T.E
        Points_Trans_Wake_P[i_1,2] = pitch_trans_wake     # Grid points for the transition wake (pitch) - T.E -
        # It is costant everywhere (right now) because V_tang does not take into
        # account of the induced velocity (due the fact that we don't know it yet)

        delta_trans_wake = (-4 * Var.Rad_P - x_trans_wake)/(Var.N_P_L)
        # The transition wake goes four radii downstream
        # Loop used to divide the transition wake in N_P_L parts (N_P_L+1 points)
        for j in range(Var.N_P_L):
            i_2 = (i_1) + j+1

            # Grid points for the transition wake
            Points_Trans_Wake_P[i_2,0] = x_trans_wake + (j+1) * delta_trans_wake # (x)
            Points_Trans_Wake_P[i_2,1] = r_trans_wake                            # (radius)
            Points_Trans_Wake_P[i_2,2] = pitch_trans_wake                        # (pitch)

    with open("output/Propeller_Points_Trans_Wake_Old.txt", "w") as file:
        file.write(f"{'Point':<8}{'x':<12}{'r':<20}{'p':<20}\n")
        for i in range(Var.Msp+1):
            i_1 = i+i*(Var.N_P_L)
            file.write(f"{i:<5}{Points_Trans_Wake_P[i_1,0]:13.9f}{Points_Trans_Wake_P[i_1,1]:13.9f}"
                       f"{Points_Trans_Wake_P[i_1,2]:13.9f}\n")

            for j in range(Var.N_P_L):
                i_2 = (i_1) + j+1
                file.write(f"{i:<5}{Points_Trans_Wake_P[i_2,0]:13.9f}{Points_Trans_Wake_P[i_2,1]:13.9f}"
                           f"{Points_Trans_Wake_P[i_2,2]:13.9f}\n")

    with open("output/Propeller_Points_Trans_Wake.txt", "w") as file:
        file.write(f"{'Point':<8}{'x':<12}{'r':<20}{'p':<20}\n")
        for i in range(Var.Msp+1):
            i_1 = i+i*(Var.N_P_L)
            file.write(f"{i:<5}{Points_Trans_Wake_P[i_1,0]:13.9f}{Points_Trans_Wake_P[i_1 ,1]:13.9f}"
                       f"{Points_Trans_Wake_P[i_1,2]:13.9f}\n")

            for j in range(Var.N_P_L):
                i_2 = (i_1) + j+1
                file.write(f"{i:<5}{Points_Trans_Wake_P[i_2,0]:13.9f}{Points_Trans_Wake_P[i_2,1]:13.9f}"
                           f"{Points_Trans_Wake_P[i_2,2]:13.9f}\n")

    return(S_Distr_P, r_R_P, t_gp_P, s_gp_P, Grid_Points_P, Control_Points_P,
               N_Panel_P, N_Bound_Vortex_P, Horseshoe_P, Points_Trans_Wake_P)


(S_Distr_P, r_R_P, t_gp_P, s_gp_P, Grid_Points_P, Control_Points_P, N_Panel_P, N_Bound_Vortex_P, Horseshoe_P, Points_Trans_Wake_P
       )=Grid_Generation_Propeller()