"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This subroutine is tasked with creating and solving the system of
equations for the propeller analysis. It aligns the wake of the propeller,
ensuring that the flow dynamics are accurately represented and optimized.
"""


import numpy as np
import sources.Variables as Var
from sources.Weight_Function_Propeller_P import Weight_function_propeller
from sources.Onset_Flow_Propeller_P import Onset_Flow_Propeller
from sources.Mid_Vect_Propeller_P import Mid_Vect_Propeller
from sources.Induced_Grid_Propeller_P import Induced_Grid_Propeller
from sources.Velocity_Total_No_Onset_Propeller_P import Velocity_Total_No_Onset_Propeller
from sources.Gamma_Initialization import Gamma_It
from sources.Propeller_Pitch import pitch
from sources.Velocity_Total_Propeller_P import Velocity_Total_Propeller
from sources.Align_Wake_Propeller_P import Align_Wake_Propeller
from sources.Skin_Friction_Drag_P import Skin_Friction_Drag


def System_Equations_Propeller_P():
    Weight_P = Weight_function_propeller()
    Gamma_TE_P = Gamma_It()
    V_Onset_P = Onset_Flow_Propeller()
    V_Ind_P, V_Tral_P = Velocity_Total_No_Onset_Propeller()
    Points_Trans_Wake_P = np.loadtxt("output/Propeller_Points_Trans_Wake.txt",
                                     skiprows= 1, usecols= (1,2,3))

    # DECLARATION OF VARIABLES

    matr_T = np.zeros((Var.Msp+1, Var.Msp+1))
    matr_Q1 = np.zeros((Var.Msp+1, Var.Msp+1))
    matr_Q2 = np.zeros((Var.Msp+1, Var.Msp+1))
    matrix = np.zeros((Var.Msp+1, Var.Msp+1))
    rhsQ = np.zeros((Var.Msp+1,2))
    #Right hand side of the equation system
    rhs = np.zeros((Var.Msp+1,1))

    T_fr_P = 0.0
    # Thrust - Skin friction drag - Propeller
    Tr_P = 0.0

    R_Circ_P = np.zeros((Var.Msp))
    # Radius where the circulation is calculated at the T.E.
    R_Circ_P_R = np.zeros((Var.Msp))
    # Dimensionaless radius where the circulation is calculated at the T.E.
    pitch_0 = np.zeros((Var.Msp+1,1))

    # INITIALIZATION

    Cs_T_r = (Var.Tr_P)/Var.rho/float(Var.Z_Blade_P)
    # Required thrust for each blade without rho (We don´t use rho in the system)

    iteration = 1

    # INITIALIZATION VARIABLE GAMMA

    Gamma_TE_P_No_dim = np.zeros((Var.Msp,1))
    #Distribution of circulation at the T.E. (Dimensionaless)
    Gamma_Panel_P = np.zeros((Var.Msp*Var.Nch))
    # Distribution of circulation on the blade

    # ALIGNMENT LOOP

    for j in range (Var.Msp+1):
        rhs[j,0] = 0.0
        rhsQ[j,0] = 0.0
        rhsQ [j,1] = 0.0
        for i in range (Var.Msp+1):
            matr_T[j,i] = 0.0
            matr_Q1[j,i] = 0.0
            matr_Q2[j,i] = 0.0
            matrix[j,i] = 0.0

    # SYSTEM OF EQUATIONS - DOUBLE LOOP USED TO CALCULATE &T(Uo),&Q1(Uo),&Q2(Uo)

    # This loop creates the system of equations (m = 1,2,3... Msp - Lines of the matrix)
    for m in range(Var.Msp):
        temp_T_0  = 0.0      # Initialization of the temporary variable used to calculate &T(Uo)
        temp_Q1_0 = 0.0      # Initialization of the temporary variable used to calculate &Q1(Uo)
        temp_Q2_0 = 0.0      # Initialization of the temporary variable used to calculate &Q2(Uo)

        for n in range (Var.Nch):       # First loop used to calculate the first sum (Nch)
            npln =  (n)+(m)*Var.Nch		# Counter used to select the right panel

            n_side = 4        # Number of sides for each panel
            if n == 0:
                n_side = 3
                # If we are considering the T.E. panel, instead of removing
                # the value of the T.E. side, we skip it

            temp_T_1 = 0.0     # Initialization of the temporary variable used to calculate &T
            temp_Q_11 = 0.0    # Initialization of the temporary variable used to calculate &Q1
            temp_Q_22 = 0.0    # Initialization of the temporary variable used to calculate &Q2

            for l in range(n_side):
                # Second loop used to calculate the second sum (4)
                xxn,xyn,xzn,xln,yln,zln = Mid_Vect_Propeller(npln,l)
                # This subroutine is used to calculate the midpoint

                temp_T_1 = temp_T_1 + zln*V_Onset_P[npln,l,1] - yln*V_Onset_P[npln,l,2]
                # Temporary variable used to calculate &T

                temp_Q_11 = temp_Q_11 + xyn*yln*V_Onset_P[npln,l,0] - xyn*xln*V_Onset_P[npln,l,1]
                # Temporary variable used to calculate &Q1

                temp_Q_22 = temp_Q_22 + xzn*xln*V_Onset_P[npln,l,2] - xzn*zln*V_Onset_P[npln,l,0]
                # Temporary variable used to calculate &Q2

            temp_T_0  = temp_T_0 + Weight_P[m,n] * temp_T_1
            temp_Q1_0 = temp_Q1_0 + Weight_P[m,n] * temp_Q_11
            temp_Q2_0 = temp_Q2_0 + Weight_P[m,n] * temp_Q_22
            # Temporary variable used to calculate T(Uo) (Nch Loop)
            # Temporary variable used to calculate Q1(Uo) (Nch Loop)
            # Temporary variable used to calculate Q2(Uo) (Nch Loop)

        matr_T [m,Var.Msp] = temp_T_0   # Value of &T(Uo) in the right position in the matrix (Temporary matrix matr_T)
        rhsQ[m,0] = - temp_Q1_0         # Value of &Q1(Uo) in the right position in the matrix (Temporary matrix rhs_Q)
        rhsQ[m,1] = - temp_Q2_0         # Value of &Q2(Uo) in the right position in the matrix (Temporary matrix rhs_Q)

    # Double loop used to calculate &T(Gam),&Q1(Gam),&Q2(Gam)

    # Loop used to select the line of the equation
    # (We don´t have the loop n because we already did that in Induced_Grid_Propeller)
    for m in range(Var.Msp):
        # Loop used to select the spanwise layer that induces velocity (Columns of the matrix) - Msp SUM
        for j in range (Var.Msp):
            temp_T_Gam = 0.0          # Initialization of the temporary variable used to calculate &T(Gam)
            temp_Q1_Gam = 0.0         # Initialization of the temporary variable used to calculate &T(Gam)
            temp_Q2_Gam = 0.0         # Initialization of the temporary variable used to calculate &T(Gam)

            #Loop used to select where the point is located (chordwise) - First SUM Nch
            for n in range(Var.Nch):
                npln = n + (m)*Var.Nch        #Panel where the point is located

                temp_T_1 = 0.0         # Initialization of the temporary variable used to calculate &T
                temp_Q_11 = 0.0        # Initialization of the temporary variable used to calculate &Q1
                temp_Q_22 = 0.0		   # Initialization of the temporary variable used to calculate &Q2

                n_side = 4             # If we are considering the T.E. panel,
                if n == 0:             # instead of removing the value of the T.E. side, we skip it
                    n_side = 3

                for l in range (n_side):
                    # Loop used to select the side of the panel

                    xxn,xyn,xzn,xln,yln,zln =  Mid_Vect_Propeller(npln,l)
                    # This subroutine is used to calculate the midpoint

                    temp_T_1 = temp_T_1 + zln*V_Ind_P[j,npln,l,1] - yln*V_Ind_P[j,npln,l,2]
                    # Temporary variable used to calculate &T - Total thrust for that panel by j

                    temp_Q_11 = temp_Q_11 + xyn*yln*V_Ind_P[j,npln,l,0]- xyn*xln*V_Ind_P[j,npln,l,1]
                    # Temporary variable used to calculate &Q1 - Total torque 1 for that panel by j

                    temp_Q_22 = temp_Q_22 + xzn*xln*V_Ind_P[j,npln,l,2]- xzn*zln*V_Ind_P[j,npln,l,0]
                    # Temporary variable used to calculate &Q2 - Total torque 2 for that panel by j

                temp_T_Gam = temp_T_Gam + Weight_P[m,n] *  temp_T_1
                temp_Q1_Gam = temp_Q1_Gam + Weight_P[m,n] * temp_Q_11
                temp_Q2_Gam = temp_Q2_Gam + Weight_P[m,n] * temp_Q_22
                # Temporary variable used to calculate Q1 (Nch Loop)
                # Temporary variable used to calculate Q2 (Nch Loop)
                # Temporary variable used to calculate T (Nch Loop)

            for i in range (Var.Nch):
            # Loop used to select where the point is located (chordwise)
            # Second SUM Nch
                npli = i + (j)* Var.Nch

                temp_T_1 = 0.0
                temp_Q_11 = 0.0
                temp_Q_22 = 0.0
                # Initialization of the temporary variable used to calculate &T
                # Initialization of the temporary variable used to calculate &Q1
                # Initialization of the temporary variable used to calculate &Q2

                n_side = 4
                if i == 0:
                    n_side = 3
                    # If we are considering the T.E. panel,
                    # instead of removing the value of the T.E. side, we skip it

                for l in range(n_side):
                    xxi,xyi,xzi,xli,yli,zli = Mid_Vect_Propeller(npli,l)
                    # This subroutine is used to calculate the midpoint

                    temp_T_1 = temp_T_1 + zli*V_Ind_P[m,npli,l,1]- yli*V_Ind_P[m,npli,l,2]
                    # Temporary variable used to calculate &T
                    temp_Q_11 = temp_Q_11 + xyi*yli*V_Ind_P[m,npli,l,0] - xyi*xli*V_Ind_P[m,npli,l,1]
                    # Temporary variable used to calculate &Q1
                    temp_Q_22 = temp_Q_22 + xzi*xli*V_Ind_P[m,npli,l,2]- xzi*zli*V_Ind_P[m,npli,l,0]
                    # Temporary variable used to calculate &Q2

                temp_T_Gam = temp_T_Gam + Weight_P[j,i] * temp_T_1
                temp_Q1_Gam = temp_Q1_Gam + Weight_P[j,i] * temp_Q_11
                temp_Q2_Gam = temp_Q2_Gam + Weight_P[j,i] * temp_Q_22
                # Temporary variable used to calculate T (Nch Loop)
                # Temporary variable used to calculate Q1 (Nch Loop)
                # Temporary variable used to calculate Q2 (Nch Loop)

            matr_T[m,j] = temp_T_Gam
            matr_Q1[m,j] = temp_Q1_Gam
            matr_Q2[m,j] = temp_Q2_Gam
            # Value of &T(Gam) in the right position in the matrix (Temporary matrix matr_T)
            # Value of &Q1(Gam) in the right position in the matrix (Temporary matrix matr_T)
            # Value of &Q2(Gam) in the right position in the matrix (Temporary matrix matr_T)

    # SYSTEM OF EQUATIONS - LOOP

    V_Tot_P, V_Tot_No_Onset_P = Velocity_Total_Propeller ()
    # It is used it in order to update V_Tot_P with the new values of gamma

    # Loop for the T.E. panels (They don´t have the weight function)
    for m in range(Var.Msp):
        npl0 = (m)*Var.Nch
        n_side = 3
        temp_T_Gam = 0.0

        for l in range(n_side):
            xxm,xym,xzm,xlm,ylm,zlm = Mid_Vect_Propeller(npl0,l)
            # This subroutine is used to calculate the midpoint

            temp_T_Gam = temp_T_Gam + zlm*V_Tot_P[npl0,l,1] - ylm*V_Tot_P[npl0,l,2]
            # Temporary variable used to calculate T (Nch Loop)

        matr_T[Var.Msp,m] = temp_T_Gam
        # Value of &T (T.E.) in the right position in the matrix (Temporary matrix matr_T)

        # Loop for the other panels (They don´t have the weight function)
        for n in range(1, Var.Nch):
            npl1 = n + (m)*Var.Nch
            n_side = 4
            temp_T_2_Gam = 0.0

            # Loop for the other panels
            for l in range(n_side):
                xxm,xym,xzm,xlm,ylm,zlm =  Mid_Vect_Propeller(npl1,l)
                # This subroutine is used to calculate the midpoint
                temp_T_2_Gam = temp_T_2_Gam + zlm*V_Tot_P[npl1,l,1] - ylm*V_Tot_P[npl1,l,2]
                # Temporary variable used to calculate T (Nch Loop)
            matr_T[Var.Msp,m] = matr_T[Var.Msp,m] + Weight_P[m,n]*temp_T_2_Gam
            # Value of &T in the right position in the matrix (Temporary matrix matr_T)

    # CREATION OF THE MATRIX

    lamba_t_1 = Gamma_TE_P[Var.Msp]
    # Lagrange multiplier lambda t-1

    for i in range (Var.Msp):
        rhs[i,0] = rhsQ[i,0] - rhsQ[i,1]
        # rhs matrix

        matrix[i,Var.Msp] = matr_T[i,Var.Msp]		# System of equation (Left Matrix)
        matrix[Var.Msp,i] = matr_T[Var.Msp,i]		# System of equation (Left Matrix)

        for j in range (Var.Msp):
            matrix[i,j] = matr_Q1[i,j] - matr_Q2[i,j] + lamba_t_1*matr_T[i,j]
            # System of equation (Left Matrix)

    rhs[Var.Msp,0] = Cs_T_r + (abs(T_fr_P))/Var.rho
    # Total thrust required (Required + Skin Friction Drag Propeller)
    matrix[Var.Msp,Var.Msp] = 0.0

    # SOLVE THE SYSTEM OF EQUATIONS

    rhs = np.linalg.solve(matrix, rhs)  #it solves the system of equations
    #Computes the “exact” solution, x, of the well-determined, i.e.,
    # full rank, linear matrix equation ax = b.

    # CONVERGENS OF THE SYSTEM

    res_0 = 0.0
    # Loop used to check if the residual is below a certain small limit
    for i in range(Var.Msp+1):
        res_1 = abs(1-Gamma_TE_P[i]/rhs[i,0])

        if res_1 > res_0:
            res_0 = res_1
        Gamma_TE_P[i] = rhs[i,0]                # New values of circulation
    with open ("output/Propeller_Gamma_TE_P.txt","w") as file:
        for i in range (Var.Msp+1):
            file.write(f"{Gamma_TE_P[i]:13.9f}\n")

    while (res_0 > Var.epsi):
        V_Tot_P, V_Tot_No_Onset_P = Velocity_Total_Propeller ()
        # It is used it in order to update V_Tot_P with the new values of gamma

        # Loop for the T.E. panels (They don´t have the weight function)
        for m in range(Var.Msp):
            npl0 = (m)*Var.Nch
            n_side = 3
            temp_T_Gam = 0.0

            for l in range(n_side):
                xxm,xym,xzm,xlm,ylm,zlm = Mid_Vect_Propeller(npl0,l)
                # This subroutine is used to calculate the midpoint

                temp_T_Gam = temp_T_Gam + zlm*V_Tot_P[npl0,l,1] - ylm*V_Tot_P[npl0,l,2]
                # Temporary variable used to calculate T (Nch Loop)

            matr_T[Var.Msp,m] = temp_T_Gam
            # Value of &T (T.E.) in the right position in the matrix (Temporary matrix matr_T)

            # Loop for the other panels (They don´t have the weight function)
            for n in range(1, Var.Nch):

                npl1 = n + (m)*Var.Nch
                n_side = 4
                temp_T_2_Gam = 0.0

                # Loop for the other panels
                for l in range(n_side):
                    xxm,xym,xzm,xlm,ylm,zlm =  Mid_Vect_Propeller(npl1,l)
                    # This subroutine is used to calculate the midpoint
                    temp_T_2_Gam = temp_T_2_Gam + zlm*V_Tot_P[npl1,l,1] - ylm*V_Tot_P[npl1,l,2]
                    # Temporary variable used to calculate T (Nch Loop)
                matr_T[Var.Msp,m] = matr_T[Var.Msp,m] + Weight_P[m,n]*temp_T_2_Gam
                # Value of &T in the right position in the matrix (Temporary matrix matr_T)

        # CREATION OF THE MATRIX

        lamba_t_1 = Gamma_TE_P[Var.Msp]
        # Lagrange multiplier lambda t-1

        for i in range (Var.Msp):
            rhs[i,0] = rhsQ[i,0] - rhsQ[i,1]
            # rhs matrix

            matrix[i,Var.Msp] = matr_T[i,Var.Msp]		# System of equation (Left Matrix)
            matrix[Var.Msp,i] = matr_T[Var.Msp,i]		# System of equation (Left Matrix)

            for j in range (Var.Msp):
                matrix[i,j] = matr_Q1[i,j] - matr_Q2[i,j] + lamba_t_1*matr_T[i,j]
                # System of equation (Left Matrix)

        rhs[Var.Msp,0] = Cs_T_r + (abs(T_fr_P))/Var.rho
        # Total thrust required (Required + Skin Friction Drag Propeller)
        matrix[Var.Msp,Var.Msp] = 0.0

        # SOLVE THE SYSTEM OF EQUATIONS

        rhs = np.linalg.solve(matrix, rhs)     #it solves the system of equations
        #Computes the “exact” solution, x, of the well-determined, i.e.,
        # full rank, linear matrix equation ax = b.

        # CONVERGENS OF THE SYSTEM

        res_0 = 0.0
        # Loop used to check if the residual is below a certain small limit
        for i in range(Var.Msp+1):
            res_1 = abs (1-Gamma_TE_P[i]/ rhs[i,0])

            if res_1 > res_0:
                res_0 = res_1
            Gamma_TE_P[i] = rhs[i,0]                # New values of circulation
        with open ("output/Propeller_Gamma_TE_P.txt","w") as file:
            for i in range (Var.Msp+1):
                file.write(f"{Gamma_TE_P[i]:13.9f}\n")

    print('Iteration Propeller Number: {}'.format(iteration), 'Circulation on the propeller at the TE:')

    for i in range (Var.Msp+1):
        print (i,Gamma_TE_P[i])

    for i in range (Var.Msp):
        j = (i)*Var.Nch

        xx,xy,xz,xl,yl,zl = Mid_Vect_Propeller(j,3)
        #This subroutine is used to calculate the midpoint px,py,pz

        R_Circ_P[i] = np.sqrt(xy*xy + xz*xz)
        R_Circ_P_R[i] = R_Circ_P[i]/Var.Rad_P

        Gamma_TE_P_No_dim[i] = (Gamma_TE_P[i]*100)/(np.pi*2*Var.Rad_P*Var.V_Ship)

    with open("output/Propeller_Gamma_TE.txt","w") as file:
        file.write("     Gamma_Dim     Gamma_No_Dim    Radius\n")
        for i in range (Var.Msp):
            file.write("{:13.9f}    {:13.9f}    {:13.9f}\n.format.  {Gamma_TE_P[i]}   {Gamma_TE_P_No_dim[i]}   {R_Circ_P[i]}\n")

    with open ("output/Propeller_Gamma_TE_P.txt","w") as file:
        for i in range (Var.Msp+1):
            file.write(f"{Gamma_TE_P[i]}\n")

    with open("output/Propeller_Print_Gamma_TE.txt","w") as file:
        for i in range(Var.Msp):
            file.write(f" {Gamma_TE_P_No_dim[i]}\n")

    with open("output/Propeller_Print_Radius_TE.txt","w") as file:
        for i in range(Var.Msp):
            file.write(f"{R_Circ_P_R[i]}\n")

    # DISTRIBUTION OF CIRCULATION AT THE REST OF THE BLADE

    for i in range(Var.Msp):
        npl_TE = i * Var.Nch
        Gamma_Panel_P[npl_TE] = Gamma_TE_P[i]

        for j in range (1,Var.Nch):
            npl = j + i * Var.Nch
            Gamma_Panel_P [npl] = Gamma_Panel_P[npl_TE] * Weight_P[i,j]

    with open ("output/Propeller_Gamma_Blade.txt","w") as file:
        file.write("  Panel         Gamma\n")
        for i in range(Var.Msp*Var.Nch):
            file.write(f"   {i:3d}       {Gamma_Panel_P[i]:13.9f}\n")

    # ALIGNMENT OF THE WAKE

    pitch_0 = pitch()

    Points_Trans_Wake_P, Grid_Points_P, Control_Points_P = Align_Wake_Propeller()
    res_0 = 0.0
    # Initialization of the residual

    # Loop used to check if the residual is below a certain small limit
    for i in range(Var.Msp +1 ):
        i_1 = i+i*(Var.N_P_L)
        res_1 = abs(1 - (Points_Trans_Wake_P[i_1,2] / pitch_0[i]))

    if res_1 > res_0:
        res_0 = res_1

    while(iteration < 15): # If the the residual is greater than epsi the loop starts again
        while (res_0 > Var.epsi):
            V_Onset_P = Onset_Flow_Propeller()
            V_Grid_P = Induced_Grid_Propeller()
            V_Ind_P, V_Tral_P = Velocity_Total_No_Onset_Propeller()
            T_fr_P, Q_fr_P = Skin_Friction_Drag()

            iteration = iteration + 1

            for j in range (Var.Msp+1):
                rhs[j,0] = 0.0
                rhsQ[j,0] = 0.0
                rhsQ [j,1] = 0.0
                for i in range (Var.Msp+1):
                    matr_T[j,i] = 0.0
                    matr_Q1[j,i] = 0.0
                    matr_Q2[j,i] = 0.0
                    matrix[j,i] = 0.0

            # SYSTEM OF EQUATIONS - DOUBLE LOOP USED TO CALCULATE &T(Uo),&Q1(Uo),&Q2(Uo)

            # This loop creates the system of equations (m = 1,2,3... Msp - Lines of the matrix)
            for m in range(Var.Msp):
                temp_T_0  = 0.0     # Initialization of the temporary variable used to calculate &T(Uo)
                temp_Q1_0 = 0.0     # Initialization of the temporary variable used to calculate &Q1(Uo)
                temp_Q2_0 = 0.0     # Initialization of the temporary variable used to calculate &Q2(Uo)

                # First loop used to calculate the first sum (Nch)
                for n in range (Var.Nch):
                    npln =  (n)+(m)*Var.Nch
                    # Counter used to select the right panel

                    n_side = 4       # Number of sides for each panel
                    if n == 0:
                        n_side = 3
                        # If we are considering the T.E. panel, instead of removing
                        # the value of the T.E. side, we skip it

                    temp_T_1 = 0.0
                    temp_Q_11 = 0.0
                    temp_Q_22 = 0.0
                    # Initialization of the temporary variable used to calculate &T
                    # Initialization of the temporary variable used to calculate &Q1
                    # Initialization of the temporary variable used to calculate &Q2

                    for l in range(n_side):
                        # Second loop used to calculate the second sum (4)
                        xxn,xyn,xzn,xln,yln,zln = Mid_Vect_Propeller(npln,l)
                        # This subroutine is used to calculate the midpoint

                        temp_T_1 = temp_T_1 + zln*V_Onset_P[npln,l,1] - yln*V_Onset_P[npln,l,2]
                        # Temporary variable used to calculate &T

                        temp_Q_11 = temp_Q_11 + xyn*yln*V_Onset_P[npln,l,0] - xyn*xln*V_Onset_P[npln,l,1]
                        # Temporary variable used to calculate &Q1

                        temp_Q_22 = temp_Q_22 + xzn*xln*V_Onset_P[npln,l,2] - xzn*zln*V_Onset_P[npln,l,0]
                        # Temporary variable used to calculate &Q2

                    temp_T_0  = temp_T_0 + Weight_P[m,n] * temp_T_1
                    temp_Q1_0 = temp_Q1_0 + Weight_P[m,n] * temp_Q_11
                    temp_Q2_0 = temp_Q2_0 + Weight_P[m,n] * temp_Q_22
                    # Temporary variable used to calculate T(Uo) (Nch Loop)
                    # Temporary variable used to calculate Q1(Uo) (Nch Loop)
                    # Temporary variable used to calculate Q2(Uo) (Nch Loop)

                matr_T [m,Var.Msp] = temp_T_0     # Value of &T(Uo) in the right position in the matrix (Temporary matrix matr_T)
                rhsQ[m,0] = - temp_Q1_0           # Value of &Q1(Uo) in the right position in the matrix (Temporary matrix rhs_Q)
                rhsQ[m,1] = - temp_Q2_0           # Value of &Q2(Uo) in the right position in the matrix (Temporary matrix rhs_Q)

            # Double loop used to calculate &T(Gam),&Q1(Gam),&Q2(Gam)

            # Loop used to select the line of the equation
            # (We don´t have the loop n because we already did that in Induced_Grid_Propeller)
            for m in range(Var.Msp):
                # Loop used to select the spanwise layer that induces velocity (Columns of the matrix) - Msp SUM
                for j in range (Var.Msp):
                    temp_T_Gam = 0.0
                    temp_Q1_Gam = 0.0
                    temp_Q2_Gam = 0.0
                    # Initialization of the temporary variable used to calculate &T(Gam)
                    # Initialization of the temporary variable used to calculate &T(Gam)
                    # Initialization of the temporary variable used to calculate &T(Gam)

                    #Loop used to select where the point is located (chordwise) - First SUM Nch
                    for n in range(Var.Nch):
                        npln = n + (m)*Var.Nch
                        #Panel where the point is located

                        temp_T_1 = 0.0
                        temp_Q_11 = 0.0
                        temp_Q_22 = 0.0
                        # Initialization of the temporary variable used to calculate &T
                        # Initialization of the temporary variable used to calculate &Q1
                        # Initialization of the temporary variable used to calculate &Q2

                        n_side = 4
                        if n == 0:
                            n_side = 3
                            # If we are considering the T.E. panel,
                            # instead of removing the value of the T.E. side, we skip it

                        for l in range (n_side):
                            # Loop used to select the side of the panel

                            xxn,xyn,xzn,xln,yln,zln =  Mid_Vect_Propeller(npln,l)
                            # This subroutine is used to calculate the midpoint

                            temp_T_1 = temp_T_1 + zln*V_Ind_P[j,npln,l,1] - yln*V_Ind_P[j,npln,l,2]
                            # Temporary variable used to calculate &T -
                            # Total thrust for that panel by j

                            temp_Q_11 = temp_Q_11 + xyn*yln*V_Ind_P[j,npln,l,0]- xyn*xln*V_Ind_P[j,npln,l,1]
                            # Temporary variable used to calculate &Q1 -
                            # Total torque 1 for that panel by j

                            temp_Q_22 = temp_Q_22 + xzn*xln*V_Ind_P[j,npln,l,2]- xzn*zln*V_Ind_P[j,npln,l,0]
                            # Temporary variable used to calculate &Q2 -
                            # Total torque 2 for that panel by j

                        temp_T_Gam = temp_T_Gam + Weight_P[m,n] *  temp_T_1
                        temp_Q1_Gam = temp_Q1_Gam + Weight_P[m,n] * temp_Q_11
                        temp_Q2_Gam = temp_Q2_Gam + Weight_P[m,n] * temp_Q_22
                        # Temporary variable used to calculate Q1 (Nch Loop)
                        # Temporary variable used to calculate Q2 (Nch Loop)
                        # Temporary variable used to calculate T (Nch Loop)

                    for i in range (Var.Nch):
                    # Loop used to select where the point is located (chordwise)
                    # Second SUM Nch
                        npli = i + (j)* Var.Nch

                        temp_T_1 = 0.0
                        temp_Q_11 = 0.0
                        temp_Q_22 = 0.0
                        # Initialization of the temporary variable used to calculate &T
                        # Initialization of the temporary variable used to calculate &Q1
                        # Initialization of the temporary variable used to calculate &Q2

                        n_side = 4                # If we are considering the T.E. panel,
                        if i == 0:                # instead of removing the value of the T.E. side, we skip it
                            n_side = 3

                        for l in range(n_side):
                            xxi,xyi,xzi,xli,yli,zli = Mid_Vect_Propeller(npli,l)
                            # This subroutine is used to calculate the midpoint

                            temp_T_1 = temp_T_1 + zli*V_Ind_P[m,npli,l,1]- yli*V_Ind_P[m,npli,l,2]
                            # Temporary variable used to calculate &T
                            temp_Q_11 = temp_Q_11 + xyi*yli*V_Ind_P[m,npli,l,0] - xyi*xli*V_Ind_P[m,npli,l,1]
                            #Temporary variable used to calculate &Q1
                            temp_Q_22 = temp_Q_22 + xzi*xli*V_Ind_P[m,npli,l,2]- xzi*zli*V_Ind_P[m,npli,l,0]
                            # Temporary variable used to calculate &Q2

                        temp_T_Gam = temp_T_Gam + Weight_P[j,i] * temp_T_1
                        temp_Q1_Gam = temp_Q1_Gam + Weight_P[j,i] * temp_Q_11
                        temp_Q2_Gam = temp_Q2_Gam + Weight_P[j,i] * temp_Q_22
                        # Temporary variable used to calculate T (Nch Loop)
                        # Temporary variable used to calculate Q1 (Nch Loop)
                        #  Temporary variable used to calculate Q2 (Nch Loop)

                    matr_T[m,j] = temp_T_Gam
                    matr_Q1[m,j] = temp_Q1_Gam
                    matr_Q2[m,j] = temp_Q2_Gam
                    # Value of &T(Gam) in the right position in the matrix (Temporary matrix matr_T)
                    # Value of &Q1(Gam) in the right position in the matrix (Temporary matrix matr_T)
                    # Value of &Q2(Gam) in the right position in the matrix (Temporary matrix matr_T)

            # SYSTEM OF EQUATIONS - LOOP

            V_Tot_P, V_Tot_No_Onset_P = Velocity_Total_Propeller ()
            # It is used it in order to update V_Tot_P with the new values of gamma

            # Loop for the T.E. panels (They don´t have the weight function)
            for m in range(Var.Msp):
                npl0 = (m)*Var.Nch
                n_side = 3
                temp_T_Gam = 0.0

                for l in range(n_side):
                    xxm,xym,xzm,xlm,ylm,zlm = Mid_Vect_Propeller(npl0,l)
                    # This subroutine is used to calculate the midpoint

                    temp_T_Gam = temp_T_Gam + zlm*V_Tot_P[npl0,l,1] - ylm*V_Tot_P[npl0,l,2]
                    # Temporary variable used to calculate T (Nch Loop)

                matr_T[Var.Msp,m] = temp_T_Gam
                # Value of &T (T.E.) in the right position in the matrix (Temporary matrix matr_T)

                # Loop for the other panels (They don´t have the weight function)
                for n in range(1, Var.Nch):

                    npl1 = n + (m)*Var.Nch
                    n_side = 4
                    temp_T_2_Gam = 0.0

                    # Loop for the other panels
                    for l in range(n_side):
                        xxm,xym,xzm,xlm,ylm,zlm =  Mid_Vect_Propeller(npl1,l)
                        # This subroutine is used to calculate the midpoint
                        temp_T_2_Gam = temp_T_2_Gam + zlm*V_Tot_P[npl1,l,1] - ylm*V_Tot_P[npl1,l,2]
                        # Temporary variable used to calculate T (Nch Loop)
                    matr_T[Var.Msp,m] = matr_T[Var.Msp,m] + Weight_P[m,n]*temp_T_2_Gam
                    # Value of &T in the right position in the matrix (Temporary matrix matr_T)

            # CREATION OF THE MATRIX

            lamba_t_1 = Gamma_TE_P[Var.Msp]        # Lagrange multiplier lambda t-1

            for i in range (Var.Msp):
                rhs[i,0] = rhsQ[i,0] - rhsQ[i,1]   # rhs matrix

                matrix[i,Var.Msp] = matr_T[i,Var.Msp]		# System of equation (Left Matrix)
                matrix[Var.Msp,i] = matr_T[Var.Msp,i]		# System of equation (Left Matrix)

                for j in range (Var.Msp):
                    matrix[i,j] = matr_Q1[i,j] - matr_Q2[i,j] + lamba_t_1*matr_T[i,j]
                    # System of equation (Left Matrix)


            rhs[Var.Msp,0] = Cs_T_r + (abs(T_fr_P))/Var.rho
            # Total thrust required (Required + Skin Friction Drag Propeller)
            matrix[Var.Msp,Var.Msp] = 0.0

            # SYSTEM OF EQUATIONS - LOOP

            rhs = np.linalg.solve(matrix, rhs)          #it solves the system of equations
            #Computes the “exact” solution, x, of the well-determined, i.e.,
            # full rank, linear matrix equation ax = b.

            # CONVERGENS OF THE SYSTEM

            res_0 = 0.0
            # Loop used to check if the residual is below a certain small limit
            for i in range(Var.Msp+1):
                res_1 = abs(1-Gamma_TE_P[i]/rhs[i,0])

                if res_1 > res_0:
                    res_0 = res_1
                Gamma_TE_P[i] = rhs[i,0]                # New values of circulation
            with open ("output/Propeller_Gamma_TE_P.txt","w") as file:
                for i in range (Var.Msp+1):
                    file.write(f"{Gamma_TE_P[i]:13.9f}\n")

            while (res_0 > Var.epsi):
                V_Tot_P, V_Tot_No_Onset_P = Velocity_Total_Propeller ()
                # It is used it in order to update V_Tot_P with the new values of gamma

                # Loop for the T.E. panels (They don´t have the weight function)
                for m in range(Var.Msp):
                    npl0 = (m)*Var.Nch
                    n_side = 3
                    temp_T_Gam = 0.0

                    for l in range(n_side):
                        xxm,xym,xzm,xlm,ylm,zlm = Mid_Vect_Propeller(npl0,l)
                        # This subroutine is used to calculate the midpoint

                        temp_T_Gam = temp_T_Gam + zlm*V_Tot_P[npl0,l,1] - ylm*V_Tot_P[npl0,l,2]
                        # Temporary variable used to calculate T (Nch Loop)

                    matr_T[Var.Msp,m] = temp_T_Gam
                    # Value of &T (T.E.) in the right position in the matrix (Temporary matrix matr_T)

                    # Loop for the other panels (They don´t have the weight function)
                    for n in range(1, Var.Nch):
                        npl1 = n + (m)*Var.Nch
                        n_side = 4
                        temp_T_2_Gam = 0.0

                        # Loop for the other panels
                        for l in range(n_side):
                            xxm,xym,xzm,xlm,ylm,zlm =  Mid_Vect_Propeller(npl1,l)
                            # This subroutine is used to calculate the midpoint
                            temp_T_2_Gam = temp_T_2_Gam + zlm*V_Tot_P[npl1,l,1] - ylm*V_Tot_P[npl1,l,2]
                            # Temporary variable used to calculate T (Nch Loop)
                        matr_T[Var.Msp,m] = matr_T[Var.Msp,m] + Weight_P[m,n]*temp_T_2_Gam
                        # Value of &T in the right position in the matrix (Temporary matrix matr_T)

                # CREATION OF THE MATRIX

                lamba_t_1 = Gamma_TE_P[Var.Msp]
                # Lagrange multiplier lambda t-1

                for i in range (Var.Msp):
                    rhs[i,0] = rhsQ[i,0] - rhsQ[i,1]
                    # rhs matrix

                    matrix[i,Var.Msp] = matr_T[i,Var.Msp]		# System of equation (Left Matrix)
                    matrix[Var.Msp,i] = matr_T[Var.Msp,i]		# System of equation (Left Matrix)

                    for j in range (Var.Msp):
                        matrix[i,j] = matr_Q1[i,j] - matr_Q2[i,j] + lamba_t_1*matr_T[i,j]
                        # System of equation (Left Matrix)

                rhs[Var.Msp,0] = Cs_T_r + (abs(T_fr_P))/Var.rho
                # Total thrust required (Required + Skin Friction Drag Propeller)
                matrix[Var.Msp,Var.Msp] = 0.0

                # SOLVE THE SYSTEM OF EQUATIONS

                rhs = np.linalg.solve(matrix, rhs)     #it solves the system of equations
                #Computes the “exact” solution, x, of the well-determined, i.e.,
                # full rank, linear matrix equation ax = b.

                # CONVERGENS OF THE SYSTEM

                res_0 = 0.0
                # Loop used to check if the residual is below a certain small limit
                for i in range(Var.Msp+1):
                    res_1 = abs(1-Gamma_TE_P[i]/rhs[i,0])

                    if res_1 > res_0:
                        res_0 = res_1
                    Gamma_TE_P[i] = rhs[i,0]                # New values of circulation
                with open ("output/Propeller_Gamma_TE_P.txt","w") as file:
                    for i in range (Var.Msp+1):
                        file.write(f"{Gamma_TE_P[i]:13.9f}\n")

            print('Iteration Propeller Number: {}'.format(iteration), 'Circulation on the propeller at the TE:')

            for i in range (Var.Msp+1):
                print (i,Gamma_TE_P[i])

            for i in range (Var.Msp):
                j = (i)*Var.Nch

                xx,xy,xz,xl,yl,zl = Mid_Vect_Propeller(j,3)
                #This subroutine is used to calculate the midpoint px,py,pz

                R_Circ_P[i] = np.sqrt(xy*xy + xz*xz)
                R_Circ_P_R[i] = R_Circ_P[i]/Var.Rad_P

                Gamma_TE_P_No_dim[i] = (Gamma_TE_P[i]*100)/(np.pi*2*Var.Rad_P*Var.V_Ship)

            with open("output/Propeller_Gamma_TE.txt","w") as file:
                file.write("     Gamma_Dim     Gamma_No_Dim    Radius\n")
                for i in range (Var.Msp):
                    file.write(f"  {Gamma_TE_P[i]}   {Gamma_TE_P_No_dim[i]}   {R_Circ_P[i]}\n")

            with open ("output/Propeller_Gamma_TE_P.txt","w") as file:
                for i in range (Var.Msp+1):
                    file.write(f"{Gamma_TE_P[i]}\n")

            with open("output/Propeller_Print_Gamma_TE.txt","w") as file:
                for i in range(Var.Msp):
                    file.write(f" {Gamma_TE_P_No_dim[i]}\n")

            with open("output/Propeller_Print_Radius_TE.txt","w") as file:
                for i in range(Var.Msp):
                    file.write(f"{R_Circ_P_R[i]}\n")

            # DISTRIBUTION OF CIRCULATION AT THE REST OF THE BLADE

            for i in range(Var.Msp):
                npl_TE = i * Var.Nch
                Gamma_Panel_P[npl_TE] = Gamma_TE_P[i]

                for j in range (1,Var.Nch):
                    npl = j + i * Var.Nch
                    Gamma_Panel_P [npl] = Gamma_Panel_P[npl_TE] * Weight_P[i,j]

            with open ("output/Propeller_Gamma_Blade.txt","w") as file:
                file.write("  Panel         Gamma\n")
                for i in range(Var.Msp*Var.Nch):
                    file.write(f"   {i:3d}       {Gamma_Panel_P[i]:13.9f}\n")

            # ALIGNMENT OF THE WAKE

            pitch_0 = pitch()

            for i in range (Var.Msp+1):
                i_1 = i+i *(Var.N_P_L)
                pitch_0[i] = Points_Trans_Wake_P[i_1,2]
                # I need to save the old value of the pitch in order
                # to evaluate the residual for the pitch distribution

            Points_Trans_Wake_P, Grid_Points_P, Control_Points_P = Align_Wake_Propeller()
            res_0 = 0.0                                   # Initialization of the residual

            # Loop used to check if the residual is below a certain small limit
            for i in range(Var.Msp+1):
                i_1 = i + i*(Var.N_P_L)
                res_1 = abs(1 - (Points_Trans_Wake_P[i_1,2] / pitch_0[i]))
                if res_1 > res_0:
                    res_0 = res_1
        break

    return Gamma_TE_P_No_dim, R_Circ_P_R
