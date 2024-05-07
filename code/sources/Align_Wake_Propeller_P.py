"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine contains a subroutine designed for optimizing
propeller flow through the generation of a new grid. It calculates induced
velocities at control points on mid-chord panels to determine a new beta
angle. This is crucial for the calculation of a new pitch, necessary for
grid generation. Pitch is interpolated at control points, with both blades
and trailers sharing the same pitch.
"""


import sources.Variables as Var
import numpy as np
from sources.Weight_Function_Propeller_P import Weight_function_propeller
from sources.Induced_Grid_Propeller_P import Induced_Grid_Propeller
from sources.Panel_Induced_Velocity_Propeller_Align import Panel_Induced_Velocity_Propeller_Align
from sources.Trailing_Vortices_Propeller_P import Trailing_Vortices_Propeller


def Align_Wake_Propeller():

    Grid_Points_P = np.zeros(((Var.Msp + 1) * (Var.Nch + 1), 3))
    Control_Points_P = np.zeros(((Var.Msp) * (Var.Nch), 3))
    Radius_cp_P = np.zeros((Var.Msp + 1))
    r_R_P,X_P,Skew_P,Chord_P,Thick_P = np.loadtxt("input/grid.txt", unpack= True)
    U_0_P, U_R_P, U_T_P = np.loadtxt("input/onset.txt",unpack=True)
    Gamma_TE_P = np.loadtxt("output/Propeller_Gamma_TE_P.txt")
    Control_Points_P = np.loadtxt('output/Propeller_Control_Points.txt')
    t_gp_P = np.loadtxt("output/Propeller_t_gp.txt", skiprows = 1)
    t_cp_P = np.loadtxt("output/Propeller_t_cp.txt", skiprows = 1)
    s_gp_P = np.loadtxt("output/Propeller_s_gp.txt", skiprows = 1)
    s_cp_P = np.loadtxt("output/Propeller_s_cp.txt", skiprows = 1)
    Weight_P = Weight_function_propeller()
    N_Bound_Vortex_P = np.loadtxt("output/Propeller_N_Bound_Vortex.txt",dtype= 'int')
    N_Bound_Vortex_P = N_Bound_Vortex_P.reshape((Var.Msp+1, 1))
    data_matrix = np.loadtxt("output/Propeller_Grid_Points_geom.txt")
    Radius_gp_P = data_matrix[:, 0]
    Chord_P_gp = data_matrix[:, 1]
    Rake_P_gp = data_matrix[:, 2]
    Skew_P_gp = data_matrix[:, 3]
    data_matrix = np.loadtxt("output/Propeller_Control_Points_geom.txt")
    Chord_P_cp = data_matrix[:, 0]
    Rake_P_cp = data_matrix[:, 1]
    Skew_P_cp = data_matrix[:, 2]

    # SUBROUTINE
    r_cp = np.zeros(Var.Msp)      # Radius in the control points of the propeller where the new pitch is computed
    tan_beta = np.zeros(Var.Msp)
    beta = np.zeros (Var.Msp)
    pitch_cp = np.zeros(Var.Msp)
    pitch_gp = np.zeros(Var.Msp+1)
    sin_b = np.zeros(Var.Msp + 1)
    cos_b = np.zeros(Var.Msp + 1)
    Theta_gp_P = np.zeros(Var.Nch + 1)
    Theta_cp_P = np.zeros(Var.Nch + 1)

    mid_point = (Var.Nch//2 + 1)

    # Loop used to select the closest control points to the midchord line (Chordwise)
    for j in range (Var.Msp):
        mid_point_cp = (mid_point + j * Var.Nch) -1

        p_x_mdp = Control_Points_P[mid_point_cp,0]     # X coordinate of the chosen control point of the propeller
        p_y_mdp = Control_Points_P[mid_point_cp,1]     # Y coordinate of the chosen control point of the propeller
        p_z_mdp = Control_Points_P[mid_point_cp,2]     # Z coordinate of the chosen control point of the propeller

        r_cp[j] = np.sqrt(p_y_mdp**2 + p_z_mdp**2)
        # Radius for the chosen control point of the propeller

        # VELOCITIES IN THE CONTROL POINTS FROM THE PANELS OF THE PROPELLER

        # Initialization of the variable used to store the induced velocity
        # from the panels of the propeller (x), (y), (z)
        u_x_panels = 0
        u_y_panels = 0
        u_z_panels = 0

        # Loop used to select the spanwise level that induces velocity
        # on the control points of the propeller
        for n in range (Var.Msp):

            # Initialization of the variable used to calculate the induced
            # velocity from the panels of the propeller (x), (y), (z)
            u_x_panels_0 = 0
            u_y_panels_0 = 0
            u_z_panels_0 = 0

            # Loop used to select the panel that induces velocity
            # on the control points of the propeller
            for m in range (Var.Nch):
                npl = m + n * Var.Nch

                u_x_temp,u_y_temp,u_z_temp = Panel_Induced_Velocity_Propeller_Align(npl,0,0,p_x_mdp,p_y_mdp,p_z_mdp)
                # Induced velocity from the selected panel on
                # the chosen control point of the propeller - No bound vortex

                u_x_panels_0 = u_x_panels_0 + Weight_P[n,m] * u_x_temp
                # Temporary variable used to calculate the induced velocity
                # from the panels of the propeller (x)
                u_y_panels_0 = u_y_panels_0 + Weight_P[n,m] * u_y_temp
                # Temporary variable used to calculate the induced velocity
                # from the panels of the propeller (y)
                u_z_panels_0 = u_z_panels_0 + Weight_P[n,m] * u_z_temp
                # Temporary variable used to calculate the induced velocity
                # from the panels of the propeller (z)

            u_x_panels = u_x_panels + Gamma_TE_P[n] * u_x_panels_0        # Induced velocity from the panels of the propeller (x)
            u_y_panels = u_y_panels + Gamma_TE_P[n] * u_y_panels_0        # Induced velocity from the panels of the propeller (y)
            u_z_panels = u_z_panels + Gamma_TE_P[n] * u_z_panels_0        # Induced velocity from the panels of the propeller (z)

        # VELOCITIES IN THE CONTROL POINTS FROM THE HORSESHOE VORTEX OF THE PROPELLER

        u_x_trail = 0
        # Initialization of the variable used to calculate the induced velocity
        # from the trailing vortices of the propeller (x)
        u_y_trail = 0
        # Initialization of the variable used to calculate the induced velocity
        # from the trailing vortices of the propeller (y)
        u_z_trail = 0
        # Initialization of the variable used to calculate the induced velocity
        # from the trailing vortices of the propeller (z)

        u_x_trail_1,u_y_trail_1,u_z_trail_1 = Trailing_Vortices_Propeller(0,p_x_mdp,p_y_mdp,p_z_mdp)
        # Induced velocity from the transition wake and from
        # the semi-infinite helicoidal vortex of the propeller (First)

        for n in range (Var.Msp):             # Loop used to select the trailing vortex that induces velocity
            n_1 = n + 1                       # on the control points of the propeller
            n_2 = (n+1) * (Var.Nch+1)

            u_x_trail_2,u_y_trail_2,u_z_trail_2 = Trailing_Vortices_Propeller(n_1,p_x_mdp,p_y_mdp,p_z_mdp)
            # Induced velocity from the transition wake and from the semi-infinite
            # helicoidal vortex of the propeller (Second) selected of the propeller

            u_x_trail = u_x_trail + Gamma_TE_P[n] * (u_x_trail_1 - u_x_trail_2)
            # Induced velocity from the horseshoe vortex of the propeller (x)
            # No bound vortex
            u_y_trail = u_y_trail + Gamma_TE_P[n] * (u_y_trail_1 - u_y_trail_2)
            # Induced velocity from the horseshoe vortex of the propeller (y)
            # No bound vortex
            u_z_trail = u_z_trail + Gamma_TE_P[n] * (u_z_trail_1 - u_z_trail_2)
            # Induced velocity from the horseshoe vortex of the propeller (z)
            # No bound vortex

            u_x_trail_1 = u_x_trail_2       # For the next loop
            u_y_trail_1 = u_y_trail_2       # For the next loop
            u_z_trail_1 = u_z_trail_2       # For the next loop

        # TOTAL INDUCED VELOCITY

        u_x_tot = u_x_trail + u_x_panels    # Total induced velocity on the propeller (x)
        u_y_tot = u_y_trail + u_y_panels    # Total induced velocity on the propeller (y)
        u_z_tot = u_z_trail + u_z_panels    # Total induced velocity on the propeller (z)

        # BETA AND PITCH AT THE CONTROL POINTS

        cos_theta = p_z_mdp / r_cp[j]
        sin_theta = p_y_mdp / r_cp[j]

        U_T_P_tot = - u_y_tot * cos_theta + u_z_tot * sin_theta
        # Total tangential induced velocity in the control points of the propeller

        U_0_P_beta = np.interp(r_cp[j],r_R_P,U_0_P)          # Wake (Axial) in the control points (s)
        U_T_P_beta = np.interp(r_cp[j],r_R_P,U_T_P)          # Wake (Tangential) in the control points (s)

        tan_beta[j] = abs(-U_0_P_beta + u_x_tot)/(Var.Omega * r_cp[j] - U_T_P_tot - U_T_P_beta)		    # New tangent beta

        beta[j] = np.arctan(tan_beta[j])

        pitch_cp[j] = tan_beta[j] * 2 * np.pi * r_cp[j]

    with open("output/Propeller_Pitch_Control_Points.txt","w")as file:
        file.write(" Spanw.     Radius          Pitch/D\n")
        for j in range (Var.Msp):
            file.write(f" {j:3d}    {r_cp[j]:13.9f}   {pitch_cp[j]/(Var.Rad_P*2):13.9f}\n")

    with open ("output/Propeller_Beta.txt","w") as file:
        file.write("    Radius           Beta\n")
        for j in range (Var.Msp):
            file.write(" {:13.9f}{:3s}{:13.9f}\n".format(r_cp[j],"", np.arctan(tan_beta[j])))

    # INTERPOLATION OF THE PITCH

    if Var.Msp == 0:
        pitch_cp[Var.Msp-1] = 0

    # This loop is used to find the values of the pitch in the grid points
    # of the propeller (No tip - No Hub)
    for i in range (1,Var.Msp):
        pitch_gp[i] = np.interp(s_gp_P[i],s_cp_P,pitch_cp)

    pitch_gp[0] = (pitch_cp[1] - pitch_cp[0])/(s_cp_P[1]-s_cp_P[0])*(s_gp_P[0] - s_cp_P[0]) + pitch_cp[0] # Pitch at the hub

    pitch_gp[Var.Msp] = (pitch_cp[Var.Msp-1] - pitch_cp[Var.Msp-2])/(s_cp_P[Var.Msp-1] - s_cp_P[Var.Msp - 2] # Pitch at the tip
        )*(s_gp_P[Var.Msp] - s_cp_P[Var.Msp - 2])+ pitch_cp[Var.Msp - 2]

    with open("output/Propeller_Pitch_Grid_Points.txt","w") as file:
        file.write(" Spanw.     Pitch\n")

        for i in range (Var.Msp+1):
            file.write(f" {i:3d}{pitch_gp[i]:13.9f}\n")

    # GRID POINTS MATRIX -  CALCULATION OF BETA(S(R)),CHORD(S),SKEW(S) AND RAKE(S)

    for i in range (Var.Msp+1):

        ipl = ((i+1)*(Var.Nch+1))-(Var.Nch+1)
        # Counter used to order the Grid Points Matrix
        p_ref_gp = np.sqrt(pitch_gp[i]**2 + (2*np.pi*Radius_gp_P[i])**2)
        # Reference pitch (It has only a radial variation)
        sin_b[i] = pitch_gp[i]/p_ref_gp
        # sin(beta)
        cos_b[i] = 2*np.pi*Radius_gp_P[i]/p_ref_gp
        # cos(beta)
        for j in range (Var.Nch+1):

            npl = (j) + ipl       # Second counter to order the Grid Points Matrix
            Theta_gp_P[j] = -Skew_P_gp[i] + (t_gp_P[j] * Chord_P_gp[i] * cos_b[i]) / Radius_gp_P[i]

            # X(s,t)
            Grid_Points_P[npl, 0] = Rake_P_gp[i] + Chord_P_gp[i] * sin_b[i] * t_gp_P[j]
            # Y(s,t)
            Grid_Points_P[npl, 1] = - Radius_gp_P[i] * np.sin(Theta_gp_P[j])
            # Z(s,t)
            Grid_Points_P[npl, 2] = Radius_gp_P[i] * np.cos(Theta_gp_P[j])

    with open('output/Propeller_Grid_Points.txt', 'w') as file:
        for i in range((Var.Nch + 1) * (Var.Msp + 1)):
            file.write(f" {Grid_Points_P[i, 0]:.9f}      {Grid_Points_P[i, 1]:.9f}      {Grid_Points_P[i, 2]:.9f}\n")

     # GRID CONTROL POINTS MATRIX - CALCULATION OF BETA(S(R)),CHORD(S),SKEW(S) AND RAKE(S)

    for i in range (Var.Msp):           # Counter used to order the Grid Control Points Matrix
        ipl = ((i+1)*(Var.Nch))-(Var.Nch)

        Radius_cp_P[i] = 0.5 * (Radius_gp_P[i] + Radius_gp_P[i+1])

        p_ref_gp = np.sqrt(pitch_cp[i]**2 + (2*np.pi*Radius_cp_P[i])**2)		# Reference pitch (It has only a radial variation)

        sin_b[i] = pitch_cp[i]/p_ref_gp                     # sin(beta)
        cos_b[i] = 2*np.pi*Radius_cp_P[i]/p_ref_gp		    # cos(beta)

        for j in range (Var.Nch):							# t Loop
            npl = (j) + ipl                                 # Second counter to order the Grid Points Matrix
            Theta_cp_P[j] = - Skew_P_cp[i] + (t_cp_P[j] * Chord_P_cp[i] * cos_b[i]) / Radius_cp_P[i]

            # X(s,t)
            Control_Points_P[npl, 0] = Rake_P_cp[i] + Chord_P_cp[i] * sin_b[i] * t_cp_P[j]
            # Y(s,t)
            Control_Points_P[npl, 1] = - Radius_cp_P[i] * np.sin(Theta_cp_P[j])
            # Z(s,t)
            Control_Points_P[npl, 2] = Radius_cp_P[i] * np.cos(Theta_cp_P[j])

    with open('output/Propeller_Control_Points.txt', 'w') as file:
        for i in range((Var.Nch) * (Var.Msp)):
            file.write(f"{Control_Points_P[i, 0]:13.9f}      {Control_Points_P[i, 1]:13.9f}      {Control_Points_P[i, 2]:13.9f}\n")

    # CREATION OF THE TRANSITION WAKE (STRAIGHT LINE VORTICES)

    Points_Trans_Wake_P = np.zeros((((Var.N_P_L+1)*(Var.Msp+1)),3))

    for i in range(Var.Msp+1):
        i_1 = i+i*(Var.N_P_L)

        x_trans_wake = Grid_Points_P[N_Bound_Vortex_P[i,0],0]           # X value for the first point of the transition wake - T.E.
        y_trans_wake = Grid_Points_P[N_Bound_Vortex_P[i,0],1]           # Y value for the first point of the transition wake - T.E.
        z_trans_wake = Grid_Points_P[N_Bound_Vortex_P[i,0],2]           # Z value for the first point of the transition wake - T.E.

        pitch_trans_wake = pitch_gp[i]                           # Pitch at the T.E. (It has only a radial variation)

        r_trans_wake = np.sqrt(y_trans_wake**2 + z_trans_wake**2)        # Radius at the T.E.

        Points_Trans_Wake_P[i_1,0] = x_trans_wake      # Grid points for the transition wake (x) - T.E
        Points_Trans_Wake_P[i_1,1] = r_trans_wake      # Grid points for the transition wake (radius) - T.E
        Points_Trans_Wake_P[i_1,2] = pitch_trans_wake  # Grid points for the transition wake (pitch) - T.E

        delta_trans_wake = (-4 * Var.Rad_P - x_trans_wake)/(Var.N_P_L)		# The transition wake goes four radii downstream

        for j in range(Var.N_P_L):   # Loop used to divide the transition wake in N_P_L parts (N_P_L+1 points)
            i_2 = (i_1) + j+1

            # Grid points for the transition wake
            Points_Trans_Wake_P[i_2,0] = x_trans_wake + (j+1) * delta_trans_wake # (x)
            Points_Trans_Wake_P[i_2,1] = r_trans_wake                            # (radius)
            Points_Trans_Wake_P[i_2,2] = pitch_trans_wake                        # (pitch)

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

    return Points_Trans_Wake_P, Grid_Points_P, Control_Points_P
