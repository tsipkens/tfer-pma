import numpy as np


# == DM2ZP ================================================================#
#   Function to calculate electric mobility from a vector of mobility diameter.
#   Author:  Timothy Sipkens, 2019-10-22
#
# Inputs:
#   d           Particle mobility diameter
#   z           Integer charge state
#   T           System temperature  (Optional)
#   P           System pressure     (Optional)
#
# Outputs:
#   Zp          Electromobility
#   B           Mechanical mobility
#
# Notes:
# 2 Some of the code is adapted from Buckley et al. (2017) and Olfert
#   laboratory.
# -------------------------------------------------------------------------#
def dm2zp(d, z, T=0.0, p=0.0):
    e = 1.6022e-19  # define electron charge [C]

    if T == 0.0 or p == 0.0:
        mu = 1.82e-5  # gas viscosity [Pa*s]
        B = Cc(d) / (3 * np.pi * mu * d)  # mechanical mobility
    else:
        S = 110.4  # temperature [K]
        T_0 = 296.15  # reference temperature [K]
        vis_23 = 1.83245e-5  # reference viscosity [kg/(m*s)]
        mu = vis_23 * ((T / T_0) ** 1.5) * ((T_0 + S) / (T + S))  # gas viscosity
        # Kim et al. (2005), ISO 15900, Eqn 3
        B = Cc(d, T, p) / (3 * np.pi * mu * d)  # mechanical mobility

    Zp = B * e * z  # electromobility

    return B, Zp


# == CC ===================================================================#
#   Function to evaluate Cunningham slip correction factor.
#   Author:  Timothy Sipkens, 2019-10-22
#
# Inputs:
#   d           Particle mobility diameter
#   T           System temperature  (Optional)
#   P           System pressure     (Optional)
#
# Outputs:
#   Zp          Electromobility
#   B           Mechanical mobility
# -------------------------------------------------------------------------#
def Cc(d, T=0.0, p=0.0):
    if T == 0.0 or p == 0.0:  # if P and T are not specified, use Buckley/Davies
        mfp = 66.5e-9  # mean free path

        # for air, from Davies (1945)
        A1 = 1.257
        A2 = 0.4
        A3 = 0.55

    else:  # from Olfert laboratory / Kim et al.
        S = 110.4  # temperature [K]
        mfp_0 = 6.730e-8  # mean free path of gas molecules in air [m]
        T_0 = 296.15  # reference temperature [K]
        p_0 = 101325  # reference pressure, [Pa] (760 mmHg to Pa)

        p = p * p_0

        mfp = mfp_0 * (T / T_0) ** 2 * (p_0 / p) * ((T_0 + S) / (T + S))  # mean free path
        # Kim et al. (2005) (doi:10.6028/jres.110.005), ISO 15900 Eqn 4

        A1 = 1.165
        A2 = 0.483
        A3 = 0.997 / 2

    Kn = (2 * mfp) / d  # Knudsen number
    Cc = 1 + Kn * (A1 + A2 * np.exp(-(2 * A3) / Kn))  # Cunningham slip correction factor

    return Cc


# == MP2ZP ================================================================#
#   Calculate electric mobility from a vector of particle mass.
#   Author:   Timothy Sipkens, 2019-10-22
#
# Inputs:
#   m           Particle mass
#   z           Integer charge state
#   T           System temperature
#   P           System pressure
#
# Outputs:
#   Zp          Electromobility
#   B           Mechanical mobility
#   d           Mobility diameter (simplied by mass-mobility relation)
#
# Note:
#   Uses mass-mobility relationship to first convert to a mobility
#   diameter and then estimates the mobility using dm2zp.
# -------------------------------------------------------------------------#
def mp2zp(m, z, T=0.0, P=0.0, prop={}):
    # -- Parse inputs ---------------------------------------------------------#
    if not bool(prop) or hasattr(prop, 'm0') or hasattr(prop, 'Dm'):
        sys.exit('Please specify the mass-mobility relation parameters in prop.')
        # if isempty or missing elements, produce error
    # -------------------------------------------------------------------------#

    d = 1e-9 * (m / prop['m0']) ** (1 / prop['Dm'])
    # use mass-mobility relationship to get mobility diameter

    # -- Use mobility diameter to get particle electro and mechanical mobl. ---#
    if T == 0.0 or P == 0.0:
        B, Zp = dm2zp(d, z)
    else:
        B, Zp = dm2zp(d, z, T, P)

    return B, Zp, d
