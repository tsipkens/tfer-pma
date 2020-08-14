
# Python tools for PMA transfer function evaluation (py-tfer-pma)

[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

The attached python functions evaluate the transfer function of particle mass analyzers (PMAs), including the centrifugal particle mass analyzer (CPMA) and aerosol particle mass analyzer (APM). This is primarily done using a novel set of expressions derived from particle tracking methods, information for which is given in an associated paper [(Sipkens, Olfert, and Rogak, 2020a)][ast20]. Further information on the different methods in this program is given as header information for each function.

This repository is intended to mirror the [mat-tfer-pma](https://github.com/tsipkens/mat-tfer-pma) repository. More information is available in the [README](https://github.com/tsipkens/mat-tfer-pma/blob/master/README.md) in that repository. 

### Dependencies

This program require *numpy*; *scipy*; and generally relies on *matplotlib* for plotting. Full function requries that these packages be installed. 

### A simple demonstration
To start, import the numpy library and, since we will be plotting data, the matplotlib.pyplot module,

```Python
import matplotlib.pyplot as plt
import numpy as np
```

Also, import the tfer_pma module, which contains the relevant functions to evaluate the transfer function using the analytical methods from [Sipkens, Olfert, and Rogak (2020a)][ast20]:

```Python
import tfer_pma # import PMA transfer function methods
```

Now define some fundamental properties, including the mass setpoint, `m_star`; the masses at which to evaluate the transfer function, `m`; the mobility diameter of the particles, `d` (note, using `d = None` will result in using the mass-mobility relation, using the values in the `prop` dictionary defined below); and the integer charge state at which to evaluate the transfer function, `z`:

```Python
# define input variables
m_star = 0.01e-18

# numpy array spanning 80 to 120% m_star
m = np.arange(0.8,1.2,0.001) * m_star
d = None # mobility diameter
z = 1. # integer charge state
```

Next, generate a dictionary that contains the properties of the particle mass analyzer, such as its geometry dimensions. Here, we also modify the default mass-mobility parameters to be used in the remainder of the program:  

```Python
prop = tfer_pma.prop_pma() # get default PMA properties

# Modify some of the properties, 
# in this case for the mass-mobility relation.
rho_eff = 900; # effective density
prop['rho0'] = rho_eff * np.pi / 6; # only used to find Rm
prop['Dm'] = 3
```

The default parameters here correspond to a centrifugal particle mass analyzer (CPMA), where the electrodes rotate at different speed. We could return an aerosol particle mass analyzer (APM) by setting the ratio of electrode speeds to unity, using:

```Python
prop['omega_hat'] = 1; # ratio of angular speeds (CPMA < 1 vs APM = 1)
```

Now we generate a setpoint dictionary. This quantity is crucial in this program, taking the `prop` dictionary generated above and two name-value pair arguments that specify the setpoint for the PMA. For example, using the mass setpoint `m_star` above and a resolution (as defined by [Reavell, Symonds, and Rushton (2011)][reavell]) of 10, we can compute the other relevant parameters to describe the PMA setpoint using:

```Python
sp,_ = tfer_pma.get_setpoint(prop, 'm_star', m_star, 'Rm', 10)
```

The output dictionary will also contain information like the voltage, `sp['V']`; angular speeds of the inner and outer electrodes, `sp['omega1']` and `sp['omega2']`, respectively; among other relevant properties. One can also use other pairings, such as specifying the voltage and centerline rotation speed: 

```Python
sp,_ = tfer_pma.get_setpoint(prop, 'V', 24.44, 'omega', 2543.9)
```

This should give a similar setpoint to the preceding statement, but specifies the setpoint in a different way. It is worth noting that most combinations of two of these parameters will be sufficient to specify to setpoint, with the exception of combining the rotational speed or voltage with a resolution, which will result in an error. 

Finally, let's evaluate the transfer function for some of the cases considered in [Sipkens, Olfert, and Rogak (2020a)][ast20]. First, consider **Case 1S**, where the fluid velocity profile is approximated using a 1st-order Taylor series expansion about the equilibrium radius. To do so: 

```Python
Lambda_1S,_ = tfer_pma.tfer_1S(sp, m, d, z, prop)
```

Here, `Lambda_1S` will be a numpy array of the same length as `m`. The other expressions from [Sipkens, Olfert, and Rogak (2020a)][ast20] can be realized using different methods from tfer_pma module, generally adopting intuitive names corresponding to the case codes from that work. For example, for **Case 1C**: 

```Python
Lambda_1C,_ = tfer_pma.tfer_1C(sp, m, d, z, prop)
```

Adding diffusion to this scenario can be done by adding `_diff` to the end of the method name above: 

```Python
Lambda_1C_diff,_ = tfer_pma.tfer_1C_diff(sp, m, d, z, prop)
```

To finish, plot the evaluate transfer functions,, with the result resembling some of the plots in [Sipkens, Olfert, and Rogak (2020a)][ast20]:

```Python
# plot the various transfer functions 
plt.plot(m, Lambda_1S)
plt.plot(m, Lambda_1C)
plt.plot(m, Lambda_1C_diff)
plt.show()
```

This sample code is given in the main.py script that is provided with this program. 

----------------------------------------------------------------------

#### License

This software is licensed under an MIT license (see the corresponding file or details).

#### Contact

This program was written by Timothy A. Sipkens ([tsipkens@mail.ubc.ca](mailto:tsipkens@mail.ubc.ca)) while at the University of British Columbia. This program was written in close association with Steven Rogak (University of British Columbia) and Jason Olfert (University of Alberta).

#### How to cite

This code should be cited by:

1. citing the associated journal article describing the particle tracking methods used in this program [(Sipkens, Olfert, and Rogak, 2020a)][ast20], and

2. optionally citing the code directly (making reference to the GitHub repository at https://github.com/tsipkens/py-tfer-pma).

#### References

[Reavell, K., J. P. R. Symonds, and M. G. Rushton. 2011. Simplified approximations to centrifugal particle mass analyser performance. Poster presented at the European Aerosol Conference, Manchester, UK, September 4.][reavell]

[Sipkens, T. A., J. S. Olfert, and S. N. Rogak. 2020a. New approaches to calculate the transfer function of particle mass analyzers. *Aerosol Sci. Technol.* 54:1, 111-127. DOI: 10.1080/02786826.2019.1680794.][ast20]

[ast20]: https://doi.org/10.1080/02786826.2019.1680794

[reavell]: https://www.researchgate.net/publication/267448365_Simplified_Approximations_to_Centrifugal_Particle_Mass_Analyser_Performance
