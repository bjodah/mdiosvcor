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
from sympy_helpers import get_sympified, get_unitless


import IAPWS95_constants as const

# Convention in this document:
# _x denotes a numerical variable e.g. 3.2353
# x_ denotes a symbol 'x'

# Variables
P_	= Symbol('P')   # Pressure (Intensive state variable)
T_	= Symbol('T')   # Temperature (Intensive state variable)
#rho_	= Symbol('rho') # Denisty (Intensive state variable)
rho_     = symbols('rho', cls=Function)
rho_     = rho_(P_,T_)

# Standard state,
val_P0    = const.standard_state_variables['P0']
val_Tdash = const.standard_state_variables['Tdash']

# Reference constants from Section 6.1 p. 428 (p. 42 in PDF)
density_units = units.kg/units.meter**3

# Reduced variables

delta_ = Symbol('delta')
tau_   = Symbol('tau')

# func_delta = symbols('delta', cls = Function)
# func_tau   = symbols('tau', cls = Function)


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
phi0_ = ln(delta_)+_n0[1]+_n0[2]*tau_+_n0[3]*ln(tau_)
for i in range(4,9):
    phi0_ += _n0[i]*ln(1-exp(-_gamma0[i]*tau_))

Psi_={}; Theta_={}; Delta_={}
for i in range(55,57):
    # Below Eq. 6.6 on p. 429 (p. 43 in PDF)
    Psi_[i]   = exp(-_C[i]*(delta_ - 1)**2 - _D[i]*(tau_ - 1)**2)
    Theta_[i] = (1 - tau_)+_A[i]*((delta_ - 1)**2)**(1/(2*_beta[i]))
    Delta_[i] = Theta_[i]**2 + _B[i]*((delta_ - 1)**2)**_a[i]


# Eq 6.6 p. 429
# Residual part
phi__r_ = 0
for i in range(1,8):
    phi__r_ += _n[i]*delta_**_d[i]*tau_**_t[i]
for i in range(8,52):
    phi__r_ += _n[i]*delta_**_d[i]*tau_**_t[i]*exp(-delta_**_c[i])
for i in range(52,55):
    phi__r_ += _n[i]*delta_**_d[i]*tau_**_t[i]*exp(-_alpha[i]*(delta_ - _epsilon[i])**2 - _beta[i]*(tau_ - _gamma[i])**2)
for i in range(55,57):
    phi__r_ += _n[i]*Delta_[i]**_b[i]*delta_*Psi_[i]

dphi__rddelta_ = diff(phi__r_, delta_)

# Different relations can be used to numerically find thermodynamic
# properties.

val_R = const.ref_variables['R']
val_rho_c = const.ref_variables['rho_c']
val_T_c   = const.ref_variables['T_c']


def get_explicit(relation, unitless=True):
    """
    Turns a relation in delta and tay into a
    relation in rho and T

    Since handling units in the hughe equation
    is quite heavy, unitless operation is default
    """
    const.ref_variables.return_unitless = unitless

    expl_delta_ = rho_/val_rho_c # Below eq 6.4 on p. 429 (p. 43 in PDF)
    expl_tau_   = val_T_c/T_     # Below eq 6.4 on p. 429 (p. 43 in PDF)
    expl_subs  = {delta_: expl_delta_, tau_: expl_tau_}

    return relation.subs(expl_subs)

# Pressure relation, Table 6.3, p. 431 (p. 45 in PDF)

pressure_relation      = P_/(rho_*val_R*T_) - 1 - delta_*dphi__rddelta_
expl_pressure_relation = get_explicit(pressure_relation)


def get_water_density(val_P=None, val_T=None, val_rho0=None, verbose = False, abstol=1e-9):
    """
    This work only uses the density of water, therefore a helper
    function is definied for accessing density at a given pressure
    and termperature
    """

    if not val_P:  val_P = val_P0
    if not val_T:  val_T = val_Tdash
    if not val_rho0: val_rho0  = sympify(1000.00) * density_units

    find_root_kwargs = {'dx0':    -1.0 * density_units,
			'xunit':  density_units,
                        'maxiter': 25,
			'xabstol': abstol,
			'verbose': verbose}

    # If P is without unit, assume Pascal:
    try:
	float(val_P/units.pascal)
    except:
	if verbose: print 'Assumeing pressure eneterd in Pascal'
	val_P *= units.pascal

    # If T is without unit, assume Kelvin:
    try:
	float(val_T/units.kelvin)
    except:
	if verbose: print 'Assumeing temperature eneterd in Kelvin'
	val_T *= units.kelvin

    # If rho0 is without unit, assume units.kg/units.meter**3:
    try:
	float(val_rho0/(units.kg/units.meter**3))
    except:
	if verbose: print 'Assuming density entered in kg / m^3'
	val_rho0 *= units.kg/units.meter**3


    return solve_relation_num(expl_pressure_relation,
				     {P_: val_P, T_: val_T},
				     rho_, val_rho0,
				     **find_root_kwargs)

def test_get_water_density(verbose=False):
    assert abs(get_unitless(get_water_density(verbose=verbose)) - 997.05) < 1e-2

def get_water_density_derivatives(diff_wrt, val_P=None, val_T=None, verbose = False, abstol=1e-9):
    if not val_P: val_P = sympify(101.3e3) * units.pascal
    if not val_T: val_T = sympify(298.15)  * units.kelvin

    # If P is without unit, assume Pascal:
    try:
        float(val_P/units.pascal)
    except:
        val_P *= units.pascal

    # If T is without unit, assume Kelvin:
    try:
        float(val_T/units.kelvin)
    except:
        val_T *= units.kelvin


    rho0=sympify(1000.00) * units.kg/units.meter**3

    find_root_kwargs = {'xunit': density_units,
			'xabstol': abstol,
			'verbose': verbose}

    return solve_relation_for_derivatives(expl_pressure_relation,
					  {P: val_P, T: val_T},
					  rho, rho0, diff_wrt,
					  **find_root_kwargs)

def test_get_water_density_derivatives(verbose=False):
    drho, drho_err = get_water_density_derivatives({T:2,P:2},
						   verbose=verbose)
    assert abs(drho[((P,1),(T,1))]-1.1913e-9)<1e-3



