#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division # so that 1/3 returns 0.333333 not 0
from collections import namedtuple
from sympy import symbols, sympify
from sympy.physics import units
from prj_helpers import ParameterStore

# Variables
# pressure and temperature (Intensive state variables)
P, T  = symbols('P, T')

# Physical constants
P0		= 101.3e3 * units.pascal
Tdash		= 298.15  * units.kelvin
eps0		= 8.854187817e-12 * units.farad / units.meter
alpha_LS	= -2.837297
M_W		= (1.00794*2+15.9994)*1e-3 * units.kilogram / units.mole
N_A		= 6.02214179e+23 / units.mole
# Page 25 in manuscript
gamma_prime_val	= 8.2e-3 * 1.602176487e-19 * 1e-18 * \
		  units.coulomb * units.meter**2


# Taylor expansion of P and T dependent function to second order
f0,dfdP,dfdT,d2fdP2,d2fdPdT,d2fdT2 = \
    symbols('f0 dfdP dfdT d2fdP2 d2fdPdT d2fdT2')
f_tay2 = f0 + dfdP*(P-P0) + dfdT*(T-Tdash) + \
       1/2*(d2fdP2*(P-P0)**2+2*d2fdPdT*(P-P0)*(T-Tdash)+d2fdT2*(T-Tdash)**2)

TayParams = namedtuple('TayParams', ['f0', 'dfdP', 'dfdT', 'd2fdP2', 'd2fdPdT', 'd2fdT2'])


# Numerical values for parameters:

# From Table 3 p. 43 in submitted manuscript

R_I_sod_params = TayParams(1.68e-10 * units.meter,
                           0        * units.meter / units.bar,
                          -4e-14    * units.meter / units.kelvin,
                           0        * units.meter / units.bar / units.bar,
                           0        * units.meter / units.bar / units.kelvin,
                           0        * units.meter / units.kelvin / units.kelvin)
R_I_sod = f_tay2.subs(R_I_sod_params._asdict())

R_I_cls_params = TayParams(2.46e-10 * units.meter,
                           0        * units.meter / units.bar,
                          -4e-14    * units.meter / units.kelvin,
                           0        * units.meter / units.bar / units.bar,
                           0        * units.meter / units.bar / units.kelvin,
                           0        * units.meter / units.kelvin / units.kelvin)
R_I_cls = f_tay2.subs(R_I_cls_params._asdict())

chi_prime_params = TayParams(7.3e-1   * units.volt,
                             0        * units.volt / units.bar,
                            -4.8e-4   * units.volt / units.kelvin,
                             0        * units.volt / units.bar / units.bar,
                             0        * units.volt / units.bar / units.kelvin,
                             0        * units.volt / units.kelvin / units.kelvin)
chi_prime = f_tay2.subs(chi_prime_params._asdict())

chi_tilde_minus_prime_params = TayParams(-1.1e-10   * units.volt * units.meter,
                             0        * units.volt * units.meter / units.bar,
                            -6.72e-14 * units.volt * units.meter / units.kelvin,
                             0        * units.volt * units.meter / units.bar / units.bar,
                             0        * units.volt * units.meter / units.bar / units.kelvin,
                             0        * units.volt * units.meter / units.kelvin / units.kelvin)
chi_tilde_minus_prime = f_tay2.subs(chi_tilde_minus_prime_params._asdict())




# Ion specific parameters
q_I_val = {'sod': sympify(1),
	   'cls': sympify(-1)}
R_I_val = {'sod': R_I_sod, # Ionic radius for each ion
	   'cls': R_I_cls} # mean value of Goldschmidt radius and
