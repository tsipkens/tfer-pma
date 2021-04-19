
import numpy as np

import re
import yaml

# == PROP_PMA =============================================================#
#   Generates the prop struct used to summarize CPMA parameters.
#   Author:   Timothy Sipkens, 2019-10-22
#
# Input:
#   opt         Options string specifying parameter set
#                   (Optional, default 'Olfert')
#
# Output:
#   prop        Properties struct for use in evaluating transfer function
# -------------------------------------------------------------------------#
def prop_pma(opts='olfert'):
    prop = {}

    with open(r'../prop/' + opts + '.yaml') as file:
        # Explicitly specify loader.
        # Allows for more variability in scientific notation.
        loader = yaml.SafeLoader
        loader.add_implicit_resolver(
            u'tag:yaml.org,2002:float',
            re.compile(u'''^(?:
             [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
            |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
            |\\.[0-9_]+(?:[eE][-+][0-9]+)?
            |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
            |[-+]?\\.(?:inf|Inf|INF)
            |\\.(?:nan|NaN|NAN))$''', re.X),
            list(u'-+0123456789.'))

        prop = yaml.load(file, Loader=loader)

    if prop == {}:
        print('Specified PMA property set is not available. Reverting to default.')
        prop = {
            'r1': 0.06,  # inner electrode radius [m]
            'r2': 0.061,  # outer electrode radius [m]
            'L': 0.2,  # length of chamber [m]
            'p': 1,  # pressure [atm]
            'T': 293,  # system temperature [K]
            'Q': 3 / 1000 / 60,  # volume flow rate (m ** 3/s) (prev: ~1 lpm)
            'omega_hat': 32 / 33  # ratio of angular speeds
        }

    # -- Parameters related to CPMA geometry ----------------------------------#
    prop['rc'] = (prop['r1'] + prop['r2']) / 2
    prop['r_hat'] = prop['r1'] / prop['r2']
    prop['del'] = (prop['r2'] - prop['r1']) / 2  # half gap width

    prop['A'] = np.pi * (prop['r2'] ** 2 - prop['r1'] ** 2)  # cross sectional area of APM
    prop['v_bar'] = prop['Q'] / prop['A']  # average flow velocity

    # -- For diffusion --------------------------------------------------------#
    kB = 1.3806488e-23  # Boltzmann's constant
    prop['D'] = lambda B: kB * prop['T'] * B  # diffusion coefficient

    # -- Default mass-mobility information -------------#
    prop['Dm'] = 3  # mass-mobility exponent
    prop['m0'] = 4.7124e-25  # mass-mobility relation density
    # Common alternate: Dm = 2.48 m0 = 2.9280e-24

    return prop
