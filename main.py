
# import functions
import matplotlib.pyplot as plt
import numpy as np

import tfer_pma


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

Lambda,G0 = tfer_pma.tfer_1S(sp, m, d, z, prop)


plt.plot(m, Lambda)
plt.show()
