"""
Date: Q4 2023 - Q1 2024
Author: Lisa Martinez
Institution: Technical University of Madrid

Description:  This subroutine is designed to open and process three specific
files containing data related to propeller characteristics, chord, and
velocities. It aims to facilitate the handling and  analysis of aerodynamic
properties through these datasets.
"""


import numpy as np
import pandas as pd
import sources.Variables as Var
from scipy.interpolate import CubicSpline


def propeller_geometry():
    # Read from exel as Data Frame
    data = pd.read_excel("input/geometry.xlsx", "geometry1", header=0)

    # Radius
    r_prop = data['Radius'].values
    # Cartesian coordinate
    x_prop = data['x'].values
    # Skew in Radians
    skew_prop = data['Skew (Rad)'].values
    # Chord
    chord_prop = data['Chord'].values
    # Thickness
    thick_prop = data['t'].values
    # Axial onset flow
    u0_prop = data['U0'].values
    # Radial onset flow
    ur_prop = data['Ur'].values
    # Tangential onset flow
    ut_prop = data['Ut'].values

    # Linear interpolation for the Radius
    ir_prop = np.linspace(r_prop[0], max(r_prop), num=Var.N_Iter)

    # Cubic spline interpolation for various properties along radial axis:
    # x coordinate, skew, chord, thickness, axial onset flow, radial onset
    # flow, tangential onset flow.
    interpolator_x = CubicSpline(r_prop, x_prop, axis=0)
    ix_prop = interpolator_x(ir_prop)

    interpolator_skew = CubicSpline(r_prop, skew_prop, axis=0)
    iskew_prop = interpolator_skew(ir_prop)

    interpolator_chord = CubicSpline(r_prop, chord_prop, axis=0)
    ichord_prop = interpolator_chord(ir_prop)

    interpolator_thick = CubicSpline(r_prop, thick_prop, axis=0)
    ithick_prop = interpolator_thick(ir_prop)

    interpolator_u0 = CubicSpline(r_prop, u0_prop, axis=0)
    iu0_prop = interpolator_u0(ir_prop)

    interpolator_ur = CubicSpline(r_prop, ur_prop, axis=0)
    iur_prop = interpolator_ur(ir_prop)

    interpolator_ut = CubicSpline(r_prop, ut_prop, axis=0)
    iut_prop = interpolator_ut(ir_prop)

    # Define datasets
    dataset = list(zip(ir_prop, ix_prop, iskew_prop, ichord_prop, ithick_prop))
    dataset1 = list(zip(iu0_prop, iur_prop, iut_prop))

    # Write dataset to 'grid.txt'
    with open('input/grid.txt', 'w') as fileID:
        format = '{:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}\n'
        for row in dataset:
            fileID.write(format.format(*row))

    # Write dataset1 to 'onset.txt'
    with open('input/onset.txt', 'w') as file:
        format = '{:10.5f} {:10.5f} {:10.5f}\n'
        for row in dataset1:
            file.write(format.format(*row))

    # Write chord_prop to 'chord.txt'
    np.savetxt('input/chord.txt', chord_prop, fmt='%10.5f')

    return ir_prop, ix_prop, iskew_prop, ichord_prop, ithick_prop


ir_prop, ix_prop, iskew_prop, ichord_prop, ithick_prop = propeller_geometry()
