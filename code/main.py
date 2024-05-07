"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description: This is the main.
"""


def main ():
    from sources.propeller_geometry import propeller_geometry

    (ir_prop, ix_prop, iskew_prop, ichord_prop, ithick_prop) = propeller_geometry()

    from sources.Grid_Generation_Propeller import Grid_Generation_Propeller
    (S_Distr_P, r_R_P, t_gp_P, s_gp_P, Grid_Points_P, Control_Points_P,
           N_Panel_P, N_Bound_Vortex_P, Horseshoe_P, Points_Trans_Wake_P
           )=Grid_Generation_Propeller()

    from sources.Weight_Function_Propeller_P import Weight_function_propeller
    Weight_P = Weight_function_propeller()

    from sources.Onset_Flow_Propeller_P import Onset_Flow_Propeller
    V_Onset_P = Onset_Flow_Propeller()

    from sources.Induced_Grid_Propeller_P import Induced_Grid_Propeller
    V_Grid_P = Induced_Grid_Propeller()

    from sources.Velocity_Total_No_Onset_Propeller_P import Velocity_Total_No_Onset_Propeller
    V_Ind_P, V_Tral_P = Velocity_Total_No_Onset_Propeller()

    from sources.System_Equations_Propeller_P import System_Equations_Propeller_P
    Gamma_TE_P_No_dim, R_Circ_P_R = System_Equations_Propeller_P()

    from sources.Advance_Ratio_P import Advance_Ratio_J
    Advance_ratio = Advance_Ratio_J()

    from sources.Skin_Friction_Drag_P import Skin_Friction_Drag
    T_fr_P, Q_fr_P = Skin_Friction_Drag()

    from sources.Efficiency_P import Efficiency
    Eff, K_T, K_Q = Efficiency()

    return Gamma_TE_P_No_dim


Gamma_TE_P_No_dim = main()
