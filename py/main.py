
# import functions
import matplotlib.pyplot as plt
import numpy as np

import tfer_pma  # import relevant functions to evaluate the PMA transfer function


# define input variables
m_star = 0.01e-18
m = np.arange(0.8, 1.2, 0.001) * m_star  # numpy array spanning 80 to 120% m_star
d = None  # mobility diameter (none uses mass-mobility relation)
z = 1.0  # integer charge state


# B,Zp,_ = tfer_pma.tfer_1S(m_star,m,d,z,prop)
prop = tfer_pma.prop_pma()  # get default PMA properties

# Modify some of the properties,
# in this case for the mass-mobility relation.
rho_eff = 900  # effective density
prop['m0'] = rho_eff * np.pi / 6 * 1e-27 # copy mass-mobility relation info (only used to find Rm)
prop['Dm'] = 3

# prop['omega_hat'] = 1; # ratio of angular speeds (CPMA < 1 vs APM = 1)

sp,_ = tfer_pma.get_setpoint(prop, 'm_star', m_star, 'Rm', 10)
# sp,_ = tfer_pma.get_setpoint(prop, 'V', 24.44, 'omega', 2543.9) # alt. phrasing


# evaluate the transfer functions
Lambda_1S, _ = tfer_pma.tfer_1S(sp, m, d, z, prop)
Lambda_1C, _ = tfer_pma.tfer_1C(sp, m, d, z, prop)
Lambda_1C_diff, _ = tfer_pma.tfer_1C_diff(sp, m, d, z, prop)
if prop['omega_hat'] == 1:
    Lambda_W1, _ = tfer_pma.tfer_W1(sp, m, d, z, prop)
    Lambda_W1_diff, _ = tfer_pma.tfer_W1_diff(sp, m, d, z, prop)

# plot the various transfer functions
plt.plot(m, Lambda_1S)
plt.plot(m, Lambda_1C)
plt.plot(m, Lambda_1C_diff)
if prop['omega_hat'] == 1:
    plt.plot(m, Lambda_W1)
    plt.plot(m, Lambda_W1_diff)
plt.show()


# generate second plot demonstrating multiple charging
m123 = np.arange(0.6,3.4,0.001) * m_star
Lambda_1C_z1, _ = tfer_pma.tfer_1C_diff(sp, m123, d, 1, prop)
Lambda_1C_z2, _ = tfer_pma.tfer_1C_diff(sp, m123, d, 2, prop)
Lambda_1C_z3, _ = tfer_pma.tfer_1C_diff(sp, m123, d, 3, prop)

plt.plot(m123, Lambda_1C_z1)
plt.plot(m123, Lambda_1C_z2)
plt.plot(m123, Lambda_1C_z3)
plt.plot(m123, Lambda_1C_z1 + Lambda_1C_z2 + Lambda_1C_z3, 'k--')  # different widths stem from resolution only applying to first pea
plt.show()
