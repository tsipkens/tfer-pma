
# import functions
import matplotlib.pyplot as plt
import numpy as np

import tfer_pma # import relevant functions to evaluate the PMA transfer function


# define input variables
m_star = 0.01e-18
m = np.arange(0.8,1.2,0.001) * m_star
d = None
z = 1.


# B,Zp,_ = tfer_pma.tfer_1S(m_star,m,d,z,prop)
prop = tfer_pma.prop_pma()

rho_eff = 900; # effective density
prop['rho0'] = rho_eff * np.pi / 6; # copy mass-mobility relation info (only used to find Rm)
prop['Dm'] = 3

sp,_ = tfer_pma.get_setpoint(prop, 'm_star', m_star, 'Rm', 10)
# sp,_ = tfer_pma.get_setpoint(prop, 'V', 24.44, 'omega', 2543.9) # alt. phrasing


# evaluate the transfer functions
Lambda_1S,_ = tfer_pma.tfer_1S(sp, m, d, z, prop)
Lambda_1C,_ = tfer_pma.tfer_1C(sp, m, d, z, prop)
Lambda_1C_diff,_ = tfer_pma.tfer_1C_diff(sp, m, d, z, prop)


# plot the various transfer functions 
plt.plot(m, Lambda_1S)
plt.plot(m, Lambda_1C)
plt.plot(m, Lambda_1C_diff)
plt.show()


