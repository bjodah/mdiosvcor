#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Solvation free energy correction terms for MD simulation of ionic solvation
# Program for calculating the correction terms applied in the publication:
#
#    Björn Dahlgren, Maria M. Reif, Philippe H. Hünenberger & Niels Hansen
#    "Calculation of derivative thermodynamic hydration and aqueous partial
#    molar properties of ions based on atomistic simulations"
#    J. Chem. Theory Comput.
#
# Latest version is available at github.com/bjodah/iosv

# Author: Björn Dahlgren
# Email (@gmail.com): bjodah
# Implemented for use in research project in the IGC group at ETH

# To the extent possible under law, Bjoern Dahlgren has waived all
# copyright and related or neighboring rights to this work.
# However, if you directly use the code in your own research, please
# cite the article referenced above.

# Any comments and/or imporvements of the code are greatly appreciated.
# Improvements are best done by making a pull request at github
# It is also the place to raise issues and filing bug reports.


# Possible future extensions (feel free to write them and make a pull
# request at github):
#   *  Implement expressions not only for LS but also for CT

from __future__ import division # 1/3 returns 0.333333 ... instead of 0

from sympy import *
from sympy import mpmath
from sympy.physics import units
# Sympy does not support PyPI package "Quantities"
from collections import namedtuple

from water_permittivity import eps
from IAPWS95_density import get_water_density

# Global variables (only LS scheme implemented)
COR_TYPES = ('B','C1','C2','D')
Y_TYPES   = ('G', 'S', 'CP', 'V', 'KT', 'AP')
# Ions (only sodium and chloride ions implemented)
IONS = ('sod','cls')

# Variables
P  = Symbol('P') # Pressure (Intensive state variable)
T  = Symbol('T') # Temperature (Intensive state variable)


# Physical constants
P0       = 101.3e3 * units.pascal
Tdash    = 298.15  * units.kelvin
#eps0     = pq.constants.vacuum_permittivity
eps0     = 8.854187817e-12 units.farad / units.meter
alpha_LS = -2.837297
M_W      = (1.00794*2+15.9994)* units.gram / units.mole
N_A      = 6.02214179e+23 / units.mole

# Taylor expansion of P and T dependent function to second order
f0,dfdP,dfdT,d2fdP2,d2fdPdT,d2fdT2 = \
    symbols('f0 dfdP dfdT d2fdP2 d2fdPdT d2fdT2')
f_tay2 = f0 + dfdP*(P-P0) + dfdT*(T-Tdash) + \
       1/2*(d2fdP2*(P-P0)**2+2*d2fdPdT*(P-P0)*(T-Tdash)+d2fdT2*(T-Tdash)**2)

# Numerical values for parameters:
TayParams = namedtuple('TayParams', ['f0', 'dfdP', 'dfdT', 'd2fdP2', 'd2fdPdT', 'd2fdT2'])

# From Table 3 p. 43 in submitted manuscript

eps_prime = eps * 66.6/get_eps(P0,Tdash)

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
q_I_val = {'sod':1,
       'cls':-1}
R_I_val = {'sod':R_I_sod, # Ionic radius for each ion
       'cls':R_I_cls} # mean value of Goldschmidt radius and

q_I, R_I = symbols('q_I, R_I')

# Dependent parameters
# Computational box length as CALCULATED Eq.37
# from density and number of water molecules
L=(N_W*M_W/N_A/rho_prime(P,T)+4*pi/3*R_I[ion](P,T)**3)**(1/3)

# Correction terms appropriate for Lattice summation
Delta_G_LS = {'B':  (8*pi*eps0)*N_A*q_I**2*(1-1/eps_prime)/L*(alpha_LS+4*pi/3*(R_I/L)**2-16*pi**2/45*(R_I/L)**5)}
	      'C1': -N_A/(6*eps0)*N_W*gamma_prime*q_I/L**3}
	      'C2': -N_A*q_I[ion]*4*pi*R_I**3/3/L**3*(chi_prime+chi_minus_prime/R_I)}
	      'D':  1/8/pi/eps0*N_A*q_I**2*(1/eps-1/eps_prime)/R_I}

Y_diff = {'G':None,
          'S':(T,1),
          'CP':(T,2),
          'V':(P,1),
          'KT':(P,2),
          'AP':(P,1,T,1)
          }

def get_Delta_Y_cor_LS(Y, P_val,T_val,N_W_val,ion,cor_type="all"):
    """
    Function which returns a correction term of type
    cor_type or a dictionary of all correction terms
    (when cor_type="all") evaluated at:

      pressure                  = P_val
      temperature               = T_val
      number of water molecules = N_W_cval
      ion type                  = ion

    Units are provided by multiplying values with a unit from
        sympy.physics.units
    If no units are provided then the pressure is assumed to be in
    Pascal and the temperature is assumed to be in Kelvin.
    """
    assert(cor_type in COR_TYPES or cor_type=="all")
    assert(ion in ions)
    assert(Y in Y_TYPES)

    if cor_type == "all":
	result = {}
	for key in Delta_G_LS.keys():
	    result[key] = get_Delta_Y_cor_LS(P_val,T_val,N_W_val,ion,cor_type=key)
	return result
    else:
        delta_Y_expr = diff(Delta_G_LS[cor_type], *Y_diff[Y])
	return delta_Y_expr.subs({P: P_val,
                                  T: T_val,
                                  N_W: N_W_val,
                                  q_I: q_I_val[ion],
                                  R_I: R_I_val[ion],
                                  })
