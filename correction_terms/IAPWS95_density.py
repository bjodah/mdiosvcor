#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Code for calculating density
# as function pressure and temperature
# according to IAPWS95 formulation.
# The code is entirely based on the publication:
#    Wagner, W., and A. Pruß. “The IAPWS Formulation 1995 for the
#    Thermodynamic Properties of Ordinary Water Substance for
#    General and Scientific Use.” Journal of Physical and Chemical
#    Reference Data 31, no. 2 (June 7, 2002): 387–535.

# Author: Björn Dahlgren
# Implemented for use in research project in the IGC group at ETH

# To the extent possible under law, Bjoern Dahlgren has waived all
# copyright and related or neighboring rights to this work.

from __future__ import division
from sympy import *
from sympy.physics import units
from root_finding import find_root, solve_relation_num, solve_relation_for_derivatives
from helpers import get_sympified

import IAPWS95_constants as const

# Standard state,
P0    =	sympify(101.325e3)*units.pascal
Tdash =	sympify(298.15)*units.kelvin


# Variables
P	= Symbol('P')   # Pressure (Intensive state variable)
T	= Symbol('T')   # Temperature (Intensive state variable)
#rho	= Symbol('rho') # Denisty (Intensive state variable)
rho     = symbols('rho', cls=Function)
rho     = rho(P,T)

# Reference constants from Section 6.1 p. 428 (p. 42 in PDF)
T_c = sympify(647.096) * units.kelvin
rho_c = sympify(322.0) * units.kg / units.meter**3
R = sympify(461.51805) * units.joule / units.kelvin / units.kg

# Reduced variables

delta = Symbol('delta')
tau   = Symbol('tau')

expl_delta = rho/rho_c # Below eq 6.4 on p. 429 (p. 43 in PDF)
expl_tau   = T_c/T     # Below eq 6.4 on p. 429 (p. 43 in PDF)
expl_subs  = {delta: expl_delta, tau: expl_tau}

func_delta = symbols('delta', cls = Function)
func_tau   = symbols('tau', cls = Function)


_n0	 = get_sympified(const._n0)
_gamma0	 = get_sympified(const._gamma0)
_n	 = get_sympified(const._n)
_c	 = get_sympified(const._c)
_d	 = get_sympified(const._d)
_t	 = get_sympified(const._t)
_alpha	 = get_sympified(const._alpha)
_beta	 = get_sympified(const._beta)
_gamma	 = get_sympified(const._gamma)
_epsilon = get_sympified(const._epsilon)
_a	 = get_sympified(const._a)
_b	 = get_sympified(const._b)
_B	 = get_sympified(const._B)
_C	 = get_sympified(const._C)
_D	 = get_sympified(const._D)
_A	 = get_sympified(const._A)


# The IAPWS-95 formulation is a fundamental equation for
# the specific Helmholtz free energy f. This equation is ex-
# pressed in dimensionless form, phi = f/(RT), and is separated
# into two parts, an ideal-gas part phi° and a residual part phi^r,
# so that:

# f(rho, T)/RT = phi = phi0 + pi__r     Eq 6.4 p. 429 (p. 43 in PDF)

# Ideal gas part, phi0
# Eq 6.5, p. 429 (p. 43 in PDF)
phi0 = ln(delta)+_n0[1]+_n0[2]*tau+_n0[3]*ln(tau)
for i in range(4,9):
    phi0 += _n0[i]*ln(1-exp(-_gamma0[i]*tau))

Psi={}; Theta={}; Delta={}
for i in range(55,57):
    # Below Eq. 6.6 on p. 429 (p. 43 in PDF)
    Psi[i]   = exp(-_C[i]*(delta - 1)**2 - _D[i]*(tau - 1)**2)
    Theta[i] = (1 - tau)+_A[i]*((delta - 1)**2)**(1/(2*_beta[i]))
    Delta[i] = Theta[i]**2 + _B[i]*((delta - 1)**2)**_a[i]


# Eq 6.6 p. 429
# Residual part
phi__r = 0
for i in range(1,8):
    phi__r += _n[i]*delta**_d[i]*tau**_t[i]
for i in range(8,52):
    phi__r += _n[i]*delta**_d[i]*tau**_t[i]*exp(-delta**_c[i])
for i in range(52,55):
    phi__r += _n[i]*delta**_d[i]*tau**_t[i]*exp(-_alpha[i]*(delta - _epsilon[i])**2 - _beta[i]*(tau - _gamma[i])**2)
for i in range(55,57):
    phi__r += _n[i]*Delta[i]**_b[i]*delta*Psi[i]

dphi__rddelta = diff(phi__r, delta)

# Pressure relation, Table 6.3, p. 431 (p. 45 in PDF)
pressure_relation      = P/(rho*R*T) - 1 - delta*dphi__rddelta
expl_pressure_relation = pressure_relation.subs(expl_subs)


# This work only uses the density of water, therefore a helper
# function is definied for accessing density at a given pressure
# and termperature



def get_water_density(P_val=None, T_val=None, verbose = False, abstol=1e-9):
    if not P_val: P_val = sympify(101.3e3) * units.pascal
    if not T_val: T_val = sympify(298.15)  * units.kelvin

    # If P is without unit, assume Pascal:
    try:
        float(P_val/units.pascal)
    except:
        P_val *= units.pascal

    # If T is without unit, assume Kelvin:
    try:
        float(T_val/units.kelvin)
    except:
        T_val *= units.kelvin


    rho0=sympify(1000.00) * units.kg/units.meter**3

    find_root_kwargs = {'dx0':    -1.0 * units.kg/units.meter**3,
			'xunit':   units.kg/units.meter**3,
			'maxiter': 25,
			'xabstol': abstol,
			'verbose': verbose}

    return solve_relation_num(expl_pressure_relation,
				     {P: P_val, T: T_val},
				     rho, rho0,
				     **find_root_kwargs)



# Entities below are currently not used. But are kept for possible
# future reuse of the code and constants
# ==================================================================
# ==================================================================



# Other constants never used explicitly
# =====================================
# Critical pressure
# P_c = 22.064e6 * units.pascal

# Triple point,
T_t	    = sympify(273.16)*units.kelvin    # Eq 2.1a p. 398
p_t	    = sympify(611.657)*units.pascal   # Eq 2.1b p. 398
rho_prime_t = sympify(999.793)*units.kg/units.meter**3
rho_bis_t   = sympify(0.00485458) * units.kg / units.meter**3



def staurated_liquid_density():
    raise NotImplemented # Not finished..
    # Saturated liquid density
    rho__prime = T*dp_sigmadT/Beta # Eq 2.3 p 398
    v = (1-T/T_c)
    b  =  [0,
           1.99274064,
           1.09965342,
          -0.510839303,
          -1.75493479,
         -45.5170352,
          -6.74694450e5]
    delta_prime = 1 + b[1]*v**(1.0/3.0) + \
                      b[2]*v**(2.0/3.0) + \
                      b[3]*v**(5.0/3.0) + \
                      b[4]*v**(16.0/3.0) + \
                      b[5]*v**(43.0/3.0) + \
                      b[6]*v**(110.0/3.0)

