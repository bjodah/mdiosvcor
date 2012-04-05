# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

from sympy import *
from sympy.physics import units
from sympy import mpmath

from __future__ import division # 1/3 returns 0.333333 ... instead of 0
#import quantities as pq #sympy does not support the python package quantitites yet

from collections import namedtuple

# <codecell>

# Variables
P  = Symbol('P') # Pressure (Intensive state variable)
T  = Symbol('T') # Temperature (Intensive state variable)

# Independent parameters (those used in the simulation)
T_val = 310.65 * units.kelvin
P_val = 5000   * units.bar
N_W   = 1024          # Number of water molecules in MD simulation
ion   = 'sod'         # What ion is used (sod: Na^+ / cls: Cl^-)

# <codecell>

# Physical constants
P0       = 101.3e3 * units.pascal
Tdash    = 298.15  * units.kelvin
#eps0     = pq.constants.vacuum_permittivity
eps0     = 8.854187817e-12 units.farad / units.meter
alpha_LS = -2.837297
M_W      = (1.00794*2+15.9994)* units.gram / units.mole
N_A      = 6.02214179e+23 / units.mole

# <codecell>

# Taylor expansion of P and T dependent function to second order
f0,dfdP,dfdT,d2fdP2,d2fdPdT,d2fdT2 = \
    symbols('f0 dfdP dfdT d2fdP2 d2fdPdT d2fdT2')
f_tay2 = f0 + dfdP*(P-P0) + dfdT*(T-Tdash) + \
       1/2*(d2fdP2*(P-P0)**2+2*d2fdPdT*(P-P0)*(T-Tdash)+d2fdT2*(T-Tdash)**2)
    
# Numerical values for parameters:
TayParams = namedtuple('TayParams', ['f0', 'dfdP', 'dfdT', 'd2fdP2', 'd2fdPdT', 'd2fdT2'])

# <codecell>

# From Table 3 p. 43
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

# <codecell>

# IAPWS 95 parametrization

# <codecell>

# Ion specific parameters
q_I = {'sod':1,
       'cls':-1}
R_I = {'sod':R_I_sod, # Ionic radius for each ion
       'cls':R_I_cls} # mean value of Goldschmidt radius and 

# <codecell>

# Dependent parameters
L = Function('L') # Computational box length as CALCULATED Eq.37
                  # from density and number of water molecules
L=(N_W*M_W/N_A/rho_prime(P,T)+4*pi/3*R_I[ion](P,T)**3)**(1/3)

eps       = Function('eps')       # Relative permitivity of REAL solvent (water)
eps_prime = Function('eps_prime') # Relative permitivity of MODEL solvent (SPC water)

# <codecell>


# <codecell>

# Correction terms appropriate for Lattice summation
Delta_G_LS_B  = (8*pi*eps0)*N_A*q_I**2*(1-1/eps_prime)/L*(alpha_LS+4*pi/3*(R_I[ion]/L)**2-16*pi**2/45*(R_I[ion]/L)**5)
Delta_G_LS_C1 = -N_A/(6*eps0)*N_W*gamma_prime*q_I/L**3
Delta_G_LS_C2 = -N_A*q_I[ion]*4*pi*R_I[ion]**3/3/L**3*(chi_prime+chi_minus_prime/R_I[ion])
Delta_G_LS_D  = 1/8/pi/eps0*N_A*q_I**2*(1/eps-1/eps_prime)/R_I[ion]
Delta_G_LS    = Delta_G_LS_B + Delta_G_LS_C1 + Delta_G_LS_C2 + Delta_G_LS_D

# <codecell>


# <codecell>


