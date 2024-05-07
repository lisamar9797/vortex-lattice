"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine computes the skin friction drag at control points
of the propeller.
"""


import numpy as np
import sources.Variables as Var
from sources.Weight_Function_Propeller_P import Weight_function_propeller
from sources.Area_Panel_P import Area_Panel
from sources.Panel_Induced_Velocity_Propeller_P import Panel_Induced_Velocity_Propeller
from sources.Trailing_Vortices_Propeller_P import Trailing_Vortices_Propeller
from sources.Biot_Savart_Propeller_P import Biot_Savart_Propeller
from sources.Mid_Vect_Propeller_P import Mid_Vect_Propeller
from sources.Normal_Vector_P import Normal_Vector


def Skin_Friction_Drag():
    I_P_Points_P = (Var.Msp*Var.Nch)
    Weight_P = Weight_function_propeller()

    # DECLARATION OF VARIABLES

    V_Tot_P = np.loadtxt("output/Propeller_Velocity_Total.txt", skiprows=2, usecols= (2,3,4))
    V_Tot_P = np.reshape(V_Tot_P, (I_P_Points_P, 4, 3))
    N_Panel_P = np.loadtxt("output/Propeller_Numeration_Panel.txt",dtype='int')
    Grid_Points_P = np.loadtxt("output/Propeller_Grid_Points.txt")
    Control_Points_P = np.loadtxt('output/Propeller_Control_Points.txt')
    Radius, beta = np.loadtxt("output/Propeller_Beta.txt",skiprows = 1, unpack = True )
    Gamma_TE_P = np.loadtxt("output/Propeller_Gamma_TE_P.txt")
    Panel, Gamma_Panel_P = np.loadtxt("output/Propeller_Gamma_Blade.txt",skiprows = 1,unpack= True)
    r_R_P,X_P,Skew_P,Chord_P,Thick_P= np.loadtxt("input/grid.txt",unpack = True)
    U_0_P, U_R_P, U_T_P = np.loadtxt("input/onset.txt",unpack = True )
    S_Distr_P,r_R_P = np.loadtxt("output/Propeller_S_Distr.txt", skiprows= 1, unpack= True)
    Points_Trans_Wake_P = np.loadtxt("output/Propeller_Points_Trans_Wake.txt",
                             skiprows= 1, usecols= (1,2,3))

    radius_cp = np.zeros((I_P_Points_P))
    s_ring = np.zeros((Var.Msp))
    vector_panel = np.zeros((3))
    tangentialDirection = np.zeros((3))
    T_Skin_F = 0.0
    Q_Skin_F = 0.0
    vector_x = np.zeros((Var.Msp))
    vector_y = np.zeros((Var.Msp))
    vector_z = np.zeros((Var.Msp))
    r_cp_a = np.zeros((Var.Msp))
    cos_theta_c = np.zeros((Var.Msp))
    sin_theta_c = np.zeros((Var.Msp))
    u_x_tot_a = np.zeros((Var.Msp))
    u_y_tot_a = np.zeros((Var.Msp))
    u_z_tot_a = np.zeros((Var.Msp))
    u_tang_skin = np.zeros((Var.Msp))
    u_rel_skin = np.zeros((Var.Msp))
    Coeff_Corr_Camber = np.zeros((Var.Msp))
    Coeff_Corr_Alpha = np.zeros((Var.Msp))
    Coeff_Corr_Thick = np.zeros((Var.Msp))
    Thick_P_skin = np.zeros((Var.Msp))
    Chord_P_skin = np.zeros((Var.Msp))
    L_Ring = np.zeros((Var.Msp))
    L_Ring_x = np.zeros((Var.Msp))
    L_Ring_y = np.zeros((Var.Msp))
    L_Ring_z = np.zeros((Var.Msp))
    C_L_Local = np.zeros((Var.Msp))
    Camber_Dimless = np.zeros((Var.Msp))
    ideal_angle_attack = np.zeros((Var.Msp))
    angle_attack = np.zeros((Var.Msp))
    beta_temp_surface = np.zeros((Var.Msp))
    pitch_cp_final_surface = np.zeros((Var.Msp))

    T_Skin_f = 0.0  # Initialization of the variable used to calculate the viscous drag
    Q_Skin_f = 0.0  # Initialization of the variable used to calculate the viscous drag

    s_tot = 0
    for j in range (Var.Msp):
        s_r = 0

        for i in range(Var.Nch):
            npl = i + (j) * Var.Nch

            x_1 = Grid_Points_P[N_Panel_P[npl,0],0]      # X value of the edge number one of the panel j
            y_1 = Grid_Points_P[N_Panel_P[npl,0],1]      # Y value of the edge number one of the panel j
            z_1 = Grid_Points_P[N_Panel_P[npl,0],2]      # Z value of the edge number one of the panel j

            x_2 = Grid_Points_P[N_Panel_P[npl,1],0]      # X value of the edge number two of the panel j
            y_2 = Grid_Points_P[N_Panel_P[npl,1],1]      # Y value of the edge number two of the panel j
            z_2 = Grid_Points_P[N_Panel_P[npl,1],2]      # Z value of the edge number two of the panel j

            x_3 = Grid_Points_P[N_Panel_P[npl,3],0]      # X value of the edge number four of the panel j
            y_3 = Grid_Points_P[N_Panel_P[npl,3],1]      # Y value of the edge number four of the panel j
            z_3 = Grid_Points_P[N_Panel_P[npl,3],2]      # Z value of the edge number four of the panel j

            x_4 = Grid_Points_P[N_Panel_P[npl,2],0]      # X value of the edge number three of the panel j
            y_4 = Grid_Points_P[N_Panel_P[npl,2],1]      # Y value of the edge number three of the panel j
            z_4 = Grid_Points_P[N_Panel_P[npl,2],2]      # Z value of the edge number three of the panel j

            s_parz = Area_Panel (x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3,x_4,y_4,z_4)
            # Area of the panel where the control point is located

            s_r = s_r + s_parz

        s_ring [j] = s_r
        s_tot = s_tot + s_r

    Ae = s_tot
    Ao = np.pi * Var.Rad_P**2

    AeAo = Ae/Ao * Var.Z_Blade_P

    # SKIN FRICTION DRAG

	# Loop used to select all the control points of the propeller
    for j in range (I_P_Points_P):

        p_x_mdp = Control_Points_P[j,0]       # X coordinate of the chosen control point of the propeller
        p_y_mdp = Control_Points_P[j,1]       # Y coordinate of the chosen control point of the propeller
        p_z_mdp = Control_Points_P[j,2]	       # Z coordinate of the chosen control point of the propeller

        x_1 = Grid_Points_P[N_Panel_P[j,0],0]     # X value of the edge number one of the panel j
        y_1 = Grid_Points_P[N_Panel_P[j,0],1]	  # Y value of the edge number one of the panel j
        z_1 = Grid_Points_P[N_Panel_P[j,0],2]     # Z value of the edge number one of the panel j

        x_2 = Grid_Points_P[N_Panel_P[j,1],0]     # X value of the edge number two of the panel j
        y_2 = Grid_Points_P[N_Panel_P[j,1],1]	  # Y value of the edge number two of the panel j
        z_2 = Grid_Points_P[N_Panel_P[j,1],2]	  # Z value of the edge number two of the panel j

        x_3 = Grid_Points_P[N_Panel_P[j,3],0]     # X value of the edge number four of the panel j
        y_3 = Grid_Points_P[N_Panel_P[j,3],1]     # Y value of the edge number four of the panel j
        z_3 = Grid_Points_P[N_Panel_P[j,3],2]     # Z value of the edge number four of the panel j

        x_4 = Grid_Points_P[N_Panel_P[j,2],0]     # X value of the edge number three of the panel j
        y_4 = Grid_Points_P[N_Panel_P[j,2],1]     # Y value of the edge number three of the panel j
        z_4 = Grid_Points_P[N_Panel_P[j,2],2]     # Z value of the edge number three of the panel j

        radius_cp[j] = np.sqrt(p_y_mdp**2 + p_z_mdp**2)  # Radius for the chosen control point of the propeller

        cos_theta_c_skin = p_z_mdp/radius_cp[j]
        sin_theta_c_skin = p_y_mdp/radius_cp[j]

        # VELOCITIES IN THE CONTROL POINTS FROM THE ONSET FLOW

        U_0_Onset = np.interp (radius_cp[j],r_R_P,U_0_P)    # Wake (Axial) in the control points (s)
        U_T_Onset = np.interp (radius_cp[j],r_R_P,U_T_P)	# Wake (Tangential) in the control points (s)
        U_R_Onset = np.interp (radius_cp[j],r_R_P,U_R_P)	# Wake (Radial) in the control points (s)

        u_x_onset = - U_0_Onset		 # Onset Flow (x)
        u_y_onset  = U_R_Onset*p_y_mdp/radius_cp[j] - U_T_Onset*p_z_mdp/radius_cp[j] + Var.Omega*p_z_mdp   # Onset Flow (y)
        u_z_onset  = U_R_Onset*p_z_mdp/radius_cp[j] + U_T_Onset*p_y_mdp/radius_cp[j] - Var.Omega*p_y_mdp   # Onset Flow (z)

        # VELOCITIES IN THE CONTROL POINTS FROM THE PANELS

        u_x_panels = 0
        # Initialization of the variable used to store the induced velocity from the panels of the propeller (x)
        u_y_panels = 0
        # Initialization of the variable used to store the induced velocity from the panels of the propeller (y)
        u_z_panels = 0
        # Initialization of the variable used to store the induced velocity from the panels of the propeller (z)

        # Loop used to select the spanwise level that induces velocity on the control points of the propeller
        for n in range (Var.Msp):
            u_x_panels_0 = 0
            #Initialization of the variable used to calculate the induced velocity from the panels of the propeller (x)
            u_y_panels_0 = 0
            # Initialization of the variable used to calculate the induced velocity from the panels of the propeller (y)
            u_z_panels_0 = 0
            # Initialization of the variable used to calculate the induced velocity from the panels of the propeller (z)

            # Loop used to select the panel that induces velocity on the control points of the propeller
            for m in range (Var.Nch):
                npl = m + (n) * Var.Nch

                u_x_temp,u_y_temp,u_z_temp = Panel_Induced_Velocity_Propeller (npl,5,0,p_x_mdp,p_y_mdp,p_z_mdp)
                # Induced velocity from the selected panel on the chosen control point of the propeller

                u_x_panels_0 = u_x_panels_0 + Weight_P[n,m] * u_x_temp
                # Temporary variable used to calculate the induced velocity from the panels of the propeller (x)
                u_y_panels_0 = u_y_panels_0 + Weight_P[n,m] * u_y_temp
                # Temporary variable used to calculate the induced velocity from the panels of the propeller (y)
                u_z_panels_0 = u_z_panels_0 + Weight_P[n,m] * u_z_temp
                # Temporary variable used to calculate the induced velocity from the panels of the propeller (z)


            u_x_panels = u_x_panels + Gamma_TE_P[n] * u_x_panels_0   # Induced velocity from the panels of the propeller (x)
            u_y_panels = u_y_panels + Gamma_TE_P[n] * u_y_panels_0   # Induced velocity from the panels of the propeller (y)
            u_z_panels = u_z_panels + Gamma_TE_P[n] * u_z_panels_0   # Induced velocity from the panels of the propeller (z)

        # VELOCITIES IN THE CONTROL POINTS FROM THE HORSESHOE VORTEX

        u_x_trail = 0
        # Initialization of the variable used to calculate the induced velocity from the trailing vortices of the propeller (x)
        u_y_trail = 0
        # Initialization of the variable used to calculate the velocity from the trailing vortices of the propeller (y)
        u_z_trail = 0
        # Initialization of the variable used to calculate the induced velocity from the trailing vortices of the propeller (z)

        x_T_E_1 = Grid_Points_P[0,0]   # First point of the first trailing vortex of the propeller (x)
        y_T_E_1 = Grid_Points_P[0,1]   # First point of the first trailing vortex of the propeller (y)
        z_T_E_1 = Grid_Points_P[0,2]   # First point of the first trailing vortex of the propeller (z)

        u_x_trail_1,u_y_trail_1,u_z_trail_1 = Trailing_Vortices_Propeller(0,p_x_mdp,p_y_mdp,p_z_mdp)
        # Induced velocity from the transition wake and from the semi-infinite helicoidal vortex of the propeller (First)

        # Loop used to select the trailing vortex that induces velocity on the control points of the propeller
        for n in range (Var.Msp):
            n_1 = n + 1
            n_2 = (n+1) * (Var.Nch+1)

            u_x_trail_2,u_y_trail_2,u_z_trail_2 = Trailing_Vortices_Propeller(n_1,p_x_mdp,p_y_mdp,p_z_mdp)
            # Induced velocity from the transition wake and from the semi-infinite helicoidal vortex of the propeller (Second)

            x_T_E_2 = Grid_Points_P[n_2,0]	# Second point of the trailing vortex of the propeller (x)
            y_T_E_2 = Grid_Points_P[n_2,1]	# Second point of the trailing vortex of the propeller (y)
            z_T_E_2 = Grid_Points_P[n_2,2]	# Second point of the trailing vortex of the propeller (z)

            U_x_s,U_y_s,U_z_s = Biot_Savart_Propeller(Var.Z_Blade_P,x_T_E_1,y_T_E_1,z_T_E_1,x_T_E_2,y_T_E_2,z_T_E_2,p_x_mdp,p_y_mdp,p_z_mdp)
            # Induced velocity from the bound vortex selected of the propeller

            u_x_trail = u_x_trail + Gamma_TE_P[n] * (u_x_trail_1 - u_x_trail_2 + U_x_s)
            # Induced velocity from the horseshoe vortex of the propeller (x)
            u_y_trail = u_y_trail + Gamma_TE_P[n] * (u_y_trail_1 - u_y_trail_2 + U_y_s)
            # Induced velocity from the horseshoe vortex of the propeller (y)
            u_z_trail = u_z_trail + Gamma_TE_P[n] * (u_z_trail_1 - u_z_trail_2 + U_z_s)
            # Induced velocity from the horseshoe vortex of the propeller (z)

            x_T_E_1 = x_T_E_2   # For the next loop
            y_T_E_1 = y_T_E_2   # For the next loop
            z_T_E_1 = z_T_E_2   # For the next loop

            u_x_trail_1 = u_x_trail_2     # For the next loop
            u_y_trail_1 = u_y_trail_2     # For the next loop
            u_z_trail_1 = u_z_trail_2     # For the next loop

        # TOTAL INDUCED VELOCITY

        u_x_tot =  u_x_onset + u_x_trail + u_x_panels  # Total induced velocity on the propeller (x)
        u_y_tot =  u_y_onset + u_y_trail + u_y_panels  # Total induced velocity on the propeller (y)
        u_z_tot =  u_z_onset + u_z_trail + u_z_panels  # Total induced velocity on the propeller (z)

        # SKIN FRICTION DRAG

        s = Area_Panel (x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3,x_4,y_4,z_4)
        # Area of the panel where the control point is located

        point_x_2,point_y_2,point_z_2,vector_xx,vector_yy,vector_zz = Mid_Vect_Propeller(j,1)
        # This subroutine is used to calculate the midpoint of the panel side number 2
        point_x_4,point_y_4,point_z_4,vector_xx,vector_yy,vector_zz = Mid_Vect_Propeller(j,3)
        # This subroutine is used to calculate the midpoint of the panel side number 2

        vector_panel[0] = point_x_4 - point_x_2    # Tangent vector to the panel (x)
        vector_panel[1] = point_y_4 - point_y_2    # Tangent vector to the panel (y)
        vector_panel[2] = point_z_4 - point_z_2    # Tangent vector to the panel (z)

        vector_panel/= np.linalg.norm(vector_panel, ord=2)     # Unit tangent vector to the panel

        tangentialDirection[0] = 0.0                 # Tangent vector in yz plane
        tangentialDirection[1] = cos_theta_c_skin    # Tangent vector in yz plane
        tangentialDirection[2]= - sin_theta_c_skin   # Tangent vector in yz plane

        tangentialDirection/= np.linalg.norm(tangentialDirection, ord = 2)   # Unit tangent vector in yz plane

        V_tang = np.dot([u_x_tot, u_y_tot, u_z_tot], vector_panel)           # Tangent velocity to the panel

        dragInPanelDirection = Var.Skin_Coeff * 0.5 * Var.rho * (abs(V_tang)*V_tang) * s  # Tangent force to the panel

        T_Skin_F = T_Skin_F + dragInPanelDirection * vector_panel[0]    # Skin friction drag (Thrust) - X force to the panel

        Q_Skin_F = Q_Skin_F + dragInPanelDirection * np.dot(vector_panel, tangentialDirection) * radius_cp[j]
            # Skin friction drag (Torque) - X force to the panel

    T_fr_P =  T_Skin_F
    Q_fr_P = - Q_Skin_F

    with open ("output/Propeller_Drag.txt", "w") as file:
        file.write("  Drag T        Drag Q          Drag KT         Drag KQ \n")
        file.write(f"{(Var.Z_Blade_P * T_fr_P):9.1f}     {(Var.Z_Blade_P * Q_fr_P):9.1f}\
        {(Var.Z_Blade_P * T_fr_P / ((Var.Omega / (2 * np.pi))**2 * Var.rho * (Var.Rad_P * 2)**4)):0.6f}\
        {Var.Z_Blade_P * Q_fr_P / ((Var.Omega / (2 * np.pi))**2 * Var.rho * (Var.Rad_P * 2)**5):0.6f}\n")

    # OPTIMIZATION PARAMETERS - OUTPUT

    mid_point = Var.Nch//2

    # Loop used to select the closest control points to the midchord line (Chordwise)
    for j in range (Var.Msp):
        mid_point_cp = (mid_point) + (j) * Var.Nch

        p_x_mdp = Control_Points_P[mid_point_cp,0]
        p_y_mdp = Control_Points_P[mid_point_cp,1]
        p_z_mdp = Control_Points_P[mid_point_cp,2]
        # X coordinate of the chosen control point of the propeller
        # Y coordinate of the chosen control point of the propeller
        # Z coordinate of the chosen control point of the propeller

        x_1 = Grid_Points_P[N_Panel_P[mid_point_cp,0],0]
        y_1 = Grid_Points_P[N_Panel_P[mid_point_cp,0],1]
        z_1 = Grid_Points_P[N_Panel_P[mid_point_cp,0],2]
        # X value of the edge number one of the panel
        # Y value of the edge number one of the panel
        # Z value of the edge number one of the panel

        x_2 = Grid_Points_P[N_Panel_P[mid_point_cp,1],0]
        y_2 = Grid_Points_P[N_Panel_P[mid_point_cp,1],1]
        z_2 = Grid_Points_P[N_Panel_P[mid_point_cp,1],2]
        # X value of the edge number two of the panel
        # Y value of the edge number two of the panel
        # Z value of the edge number two of the panel

        x_3 = Grid_Points_P[N_Panel_P[mid_point_cp,3],0]
        y_3 = Grid_Points_P[N_Panel_P[mid_point_cp,3],1]
        z_3 = Grid_Points_P[N_Panel_P[mid_point_cp,3],2]
        # X value of the edge number four of the panel
        # Y value of the edge number four of the panel
        # Z value of the edge number four of the panel

        x_4 = Grid_Points_P[N_Panel_P[mid_point_cp,2],0]
        y_4 = Grid_Points_P[N_Panel_P[mid_point_cp,2],1]
        z_4 = Grid_Points_P[N_Panel_P[mid_point_cp,2],2]
        # X value of the edge number three of the panel
        # Y value of the edge number three of the panel
        # Z value of the edge number three of the panel

        vec_x,vec_y,vec_z = Normal_Vector (x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3,x_4,y_4,z_4)
        # This subroutine calculates the normal vector for the chosen panel

        vector_x[j] = vec_x       # X component of the vector
        vector_y[j] = vec_y       # Y component of the vector
        vector_z[j] = vec_z       # Z component of the vector

        r_cp_a[j] = np.sqrt(p_y_mdp**2 + p_z_mdp**2)   # Radius for the chosen control point of the propeller

        cos_theta_c[j] = p_z_mdp/r_cp_a[j]
        sin_theta_c[j] = p_y_mdp/r_cp_a[j]

        # ONSET

        U_0_P_Onset = np.interp(r_cp_a[j],r_R_P,U_0_P)    # Wake (Axial) in the midpoints (s)
        U_T_P_Onset = np.interp(r_cp_a[j],r_R_P,U_T_P)    # Wake (Tangential) in the midpoints (s)
        U_R_P_Onset = np.interp(r_cp_a[j],r_R_P,U_R_P)	  # Wake (Radial) in the grid midpoints (s)

        u_x_onset = - U_0_P_Onset		# Onset Flow (x)
        u_y_onset = U_R_P_Onset*p_y_mdp/r_cp_a[j] - U_T_P_Onset*p_z_mdp/r_cp_a[j] + Var.Omega*p_z_mdp    # Onset Flow (y)
        u_z_onset = U_R_P_Onset*p_z_mdp/r_cp_a[j] + U_T_P_Onset*p_y_mdp/r_cp_a[j] - Var.Omega*p_y_mdp	 # Onset Flow (z)

        # VELOCITIES IN THE CONTROL POINTS FROM THE PANELS OF THE PROPELLER

        u_x_panels = 0
        # Initialization of the variable used to store the induced velocity from the panels of the propeller (x)
        u_y_panels = 0
        # Initialization of the variable used to store the induced velocity from the panels of the propeller (y)
        u_z_panels = 0
        # Initialization of the variable used to store the induced velocity from the panels of the propeller (z)

        # Loop used to select the spanwise level that induces velocity on the control points of the propeller
        for n in range (Var.Msp):
            u_x_panels_0 = 0
            # Initialization of the variable used to calculate the induced velocity from the panels of the propeller (x)
            u_y_panels_0 = 0
            # Initialization of the variable used to calculate the induced velocity from the panels of the propeller (y)
            u_z_panels_0 = 0
            # Initialization of the variable used to calculate the induced velocity from the panels of the propeller (z)

            # Loop used to select the panel that induces velocity on the control points of the propeller
            for m in range (Var.Nch):
                npl = m + (n) * Var.Nch

                u_x_temp,u_y_temp,u_z_temp = Panel_Induced_Velocity_Propeller(npl,5,0,p_x_mdp,p_y_mdp,p_z_mdp)

                # Induced velocity from the selected panel on the chosen control point of the propeller
                u_x_panels_0 = u_x_panels_0 + Weight_P[n,m] * u_x_temp
                # Temporary variable used to calculate the induced velocity from the panels of the propeller (x)
                u_y_panels_0 = u_y_panels_0 + Weight_P[n,m] * u_y_temp
                # Temporary variable used to calculate the induced velocity from the panels of the propeller (y)
                u_z_panels_0 = u_z_panels_0 + Weight_P[n,m] * u_z_temp
                # Temporary variable used to calculate the induced velocity from the panels of the propeller (z)

            u_x_panels = u_x_panels + Gamma_TE_P[n] * u_x_panels_0   # Induced velocity from the panels of the propeller (x)
            u_y_panels = u_y_panels + Gamma_TE_P[n] * u_y_panels_0   # Induced velocity from the panels of the propeller (y)
            u_z_panels = u_z_panels + Gamma_TE_P[n] * u_z_panels_0   # Induced velocity from the panels of the propeller (z)

        # VELOCITIES IN THE CONTROL POINTS FROM THE HORSESHOE VORTEX OF THE PROPELLER

        u_x_trail = 0
        # Initialization of the variable used to calculate the induced velocity from the trailing vortices of the propeller (x)
        u_y_trail = 0
        # Initialization of the variable used to calculate the induced velocity from the trailing vortices of the propeller (y)
        u_z_trail = 0
        # Initialization of the variable used to calculate the induced velocity from the trailing vortices of the propeller (z)

        x_T_E_1 = Grid_Points_P[0,0]	# First point of the first trailing vortex of the propeller (x)
        y_T_E_1 = Grid_Points_P[0,1]	# First point of the first trailing vortex of the propeller (y)
        z_T_E_1 = Grid_Points_P[0,2]	# First point of the first trailing vortex of the propeller (z)

        u_x_trail_1,u_y_trail_1,u_z_trail_1 = Trailing_Vortices_Propeller(0,p_x_mdp,p_y_mdp,p_z_mdp)
        # Induced velocity from the transition wake and from the semi-infinite helicoidal vortex of the propeller (First)

        # Loop used to select the trailing vortex that induces velocity on the control points of the propeller
        for n in range (Var.Msp):
            n_1 = n + 1
            n_2 = (n+1) * (Var.Nch+1)

            u_x_trail_2,u_y_trail_2,u_z_trail_2 = Trailing_Vortices_Propeller(n_1,p_x_mdp,p_y_mdp,p_z_mdp)
            # Induced velocity from the transition wake and from the semi-infinite helicoidal vortex of the propeller (Second)

            x_T_E_2 = Grid_Points_P[n_2,0]	# Second point of the trailing vortex of the propeller (x)
            y_T_E_2 = Grid_Points_P[n_2,1]	# Second point of the trailing vortex of the propeller (y)
            z_T_E_2 = Grid_Points_P[n_2,2]	# Second point of the trailing vortex of the propeller (z)

            U_x_s,U_y_s,U_z_s = Biot_Savart_Propeller(Var.Z_Blade_P,x_T_E_1,y_T_E_1,
                                                      z_T_E_1,x_T_E_2,y_T_E_2,z_T_E_2,p_x_mdp,p_y_mdp,p_z_mdp)
            # Induced velocity from the bound vortex selected of the propeller

            u_x_trail = u_x_trail + Gamma_TE_P[n] * (u_x_trail_1 - u_x_trail_2 + U_x_s)
            # Induced velocity from the horseshoe vortex of the propeller (x)
            u_y_trail = u_y_trail + Gamma_TE_P[n] * (u_y_trail_1 - u_y_trail_2 + U_y_s)
            # Induced velocity from the horseshoe vortex of the propeller (y)
            u_z_trail = u_z_trail + Gamma_TE_P[n] * (u_z_trail_1 - u_z_trail_2 + U_z_s)
            # Induced velocity from the horseshoe vortex of the propeller (z)

            x_T_E_1 = x_T_E_2               # For the next loop
            y_T_E_1 = y_T_E_2               # For the next loop
            z_T_E_1 = z_T_E_2               # For the next loop

            u_x_trail_1 = u_x_trail_2       # For the next loop
            u_y_trail_1 = u_y_trail_2       # For the next loop
            u_z_trail_1 = u_z_trail_2       # For the next loop

        # TOTAL INDUCED VELOCITY

        u_x_tot_a[j] = u_x_onset + u_x_trail + u_x_panels    # Total induced velocity on the propeller (x)
        u_y_tot_a[j] = u_y_onset + u_y_trail + u_y_panels    # Total induced velocity on the propeller (y)
        u_z_tot_a[j] = u_z_onset + u_z_trail + u_z_panels    # Total induced velocity on the propeller (z)

        u_tang_skin[j] = - u_y_tot_a[j] * cos_theta_c[j] + u_z_tot_a[j] * sin_theta_c[j]
        # Total tangential induced velocity in the control points of the propeller
        u_rel_skin[j] = np.sqrt((u_x_tot_a[j]**2) + (u_tang_skin[j]**2))

    # THRUST FOR EACH STRIP

    # Spanwise loop
    for m in range (Var.Msp):
        npl_TE = (m)*Var.Nch

        L_0_x = 0   # Initialization of the temporary variable used to calculate the lift of the stator (x)
        L_0_y = 0   # Initialization of the temporary variable used to calculate the lift of the stator (y)
        L_0_z = 0   # Initialization of the temporary variable used to calculate the lift of the stator (z)

        # Chordwise loop
        for n in range (Var.Nch):
            npl = n + (m)*Var.Nch

            L_00_x = 0  # Initialization of the temporary variable used to calculate the lift of the stator (x)
            L_00_y = 0  # Initialization of the temporary variable used to calculate the lift of the stator (y)
            L_00_z = 0  # Initialization of the temporary variable used to calculate the lift of the stator (z)

            # Panel loop
            for k in range (4):

                xkx,xky,xkz,xlk,ylk,zlk = Mid_Vect_Propeller(npl,k)

                L_00_x = L_00_x + zlk*V_Tot_P[npl,k,1] - ylk*V_Tot_P[npl,k,2]  # Lift generated by side k panel npl (x)
                L_00_y = L_00_y + xlk*V_Tot_P[npl,k,2] - zlk*V_Tot_P[npl,k,0]  # Lift generated by side k panel npl (y)
                L_00_z = L_00_z - xlk*V_Tot_P[npl,k,1] + ylk*V_Tot_P[npl,k,0]  # Lift generated by side k panel npl (z)

            L_0_x = L_0_x + Weight_P[m,n] * L_00_x     # Lift generated by the panel npl (x)
            L_0_y = L_0_y + Weight_P[m,n] * L_00_y     # Lift generated by the panel npl (y)
            L_0_z = L_0_z + Weight_P[m,n] * L_00_z     # Lift generated by the panel npl (z)

        xkx,xky,xkz,xlk,ylk,zlk = Mid_Vect_Propeller(npl_TE,3)

        # Lift (Ring) - x
        L_Ring_x[m] = Gamma_Panel_P[npl_TE]*L_0_x - Gamma_Panel_P[
            npl_TE]*zlk*V_Tot_P[npl_TE,3,1] + Gamma_Panel_P[npl_TE]*ylk*V_Tot_P[npl_TE,3,2]
        # Lift (Ring) - y
        L_Ring_y[m] = Gamma_Panel_P[npl_TE]*L_0_y - Gamma_Panel_P[
            npl_TE]*xlk*V_Tot_P[npl_TE,3,2] + Gamma_Panel_P[npl_TE]*zlk*V_Tot_P[npl_TE,3,0]
    	   # Lift (Ring) - z
        L_Ring_z[m] = Gamma_Panel_P[npl_TE]*L_0_z + Gamma_Panel_P[
            npl_TE]*xlk*V_Tot_P[npl_TE,3,1] - Gamma_Panel_P[npl_TE]*ylk*V_Tot_P[npl_TE,3,0]
        # Total Lift of the ring m
        L_Ring[m] = L_Ring_x[m]*vector_x[m] + L_Ring_y[m]*vector_y[m] + L_Ring_z[m]*vector_z[m]

    # CORRECTION FACTORS

    for m in range (Var.Msp):
        Chord_P_skin[m] = np.interp(r_cp_a[m],S_Distr_P,Chord_P)       # Value of the chord
        Thick_P_skin[m] = np.interp(r_cp_a[m],S_Distr_P,Thick_P)       # Value of the thickness

    Thick_0 = (Thick_P_skin[1]-Thick_P_skin[0])/(r_cp_a[1]-r_cp_a[0])*(-r_cp_a[0]) + Thick_P_skin[0]  # Pitch at the hub

    for m in range (Var.Msp):
    	Max_Skew = max(Skew_P)*180/np.pi
    	a = (3.5*AeAo) / (np.sqrt(r_cp_a[m]/Var.Rad_P * np.tan(beta[m]))) * (r_cp_a[m]/Var.Rad_P - 0.5)**2
    	b = 0.71 * np.sqrt(r_cp_a[m]/Var.Rad_P * np.tan(beta[m])) + 0.56 * (
            AeAo)**2 + (r_cp_a[m]/Var.Rad_P)*(5-Var.Z_Blade_P)/Var.Z_Blade_P + 0.46

    	Coeff_Corr_Camber[m] = a + b

    	Coeff_Corr_Thick[m] = 2*(5 + Var.Z_Blade_P)*Var.Z_Blade_P*Thick_0/(Var.Rad_P*2) * AeAo * (1-r_cp_a[m]/Var.Rad_P)**2

    	d = 1.2 * AeAo + 0.65 - 0.07*(2-np.pi*r_cp_a[m]/Var.Rad_P * np.tan(beta[m]))**3
    	e = 55/np.sqrt(np.pi*r_cp_a[m]/Var.Rad_P * np.tan(beta[m]))*AeAo*(r_cp_a[
            m]/Var.Rad_P - 0.55)**4 + 1.2*r_cp_a[m]/Var.Rad_P*(5-Var.Z_Blade_P)/Var.Z_Blade_P
    	f = 0.08 * Max_Skew * (1- 20 * abs((r_cp_a[m]/Var.Rad_P - 0.4)**3))

    	Coeff_Corr_Alpha[m]	= d + e + f

    with open("output/Propeller_Correction_Factors.txt","w") as file:
        file.write("   K_Camber   K_Thickness    K_Alpha\n")
        for m in range(Var.Msp):
            file.write(f"    {Coeff_Corr_Camber[m]:7.4f}     {Coeff_Corr_Thick[m]:7.4f}      {Coeff_Corr_Alpha[m]:7.4f}\n")

    for m in range (Var.Msp):
        C_L_Local[m] = abs(L_Ring[m])/(0.5*(u_rel_skin[m])**2*s_ring[m])      # Lift Coefficient
        Camber_Dimless[m] = (C_L_Local[m]*(1.0-Var.cny)) * 0.067 * Coeff_Corr_Camber[m] # / (1 - 0.83*Thick_P_skin[m]/Chord_P_skin[m])
        # Value of the camber
        ideal_angle_attack[m]= 1.40*(C_L_Local[m]*(1.0-Var.cny))/(180.0)*np.pi        # Ideal angle of attack - 2D
        angle_attack[m] = (C_L_Local[m]*Var.cny)/(np.pi*2.0) + ideal_angle_attack[m]  # Angle of attack -  2D
        beta_temp_surface[m] = angle_attack[m]*Coeff_Corr_Alpha[m] + beta[m] + Coeff_Corr_Thick[m]/180*np.pi
        pitch_cp_final_surface[m] = np.tan(beta_temp_surface[m]) * 2.0 * np.pi * r_cp_a[m]

    with open ("output/Propeller_Lift_Coefficient_Local_Parameters.txt","w") as file:
        file.write("  C_L_Local     Area         Beta         Angle of Attack        f/c       Ideal angle of attack      Radius\n")
        for m in range (Var.Msp):
                       file.write(f"  {C_L_Local[m]:7.4f}     {s_ring[m]:7.4f}   {beta[m]:13.9f}    {angle_attack[m]:13.9f}"
                        f"   {Camber_Dimless[m]:13.9f}      {ideal_angle_attack[m]:13.9f}     {r_cp_a[m]:13.9}\n")
    with open ("output/Propeller_Pitch_Control_Points_Angle_Add.txt","w") as file:
        file.write(" Spanw.     Radius          Pitch/D\n")
        for j in range (Var.Msp):
            file.write(f" {j:3d}    {r_cp_a[j]:13.9f}     {(pitch_cp_final_surface[j]/(Var.Rad_P*2)):13.9f}\n")

    return T_fr_P, Q_fr_P


T_fr_P,Q_fr_P = Skin_Friction_Drag()
