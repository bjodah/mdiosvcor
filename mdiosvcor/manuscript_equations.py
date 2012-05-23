#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This work is open source and is released under the
# 2-clause BSD license (see LICENSE.txt for further information)
# Copyright (c) 2011, 2012, Björn Ingvar Dahlgren

"""
Equations used when calculating the correction terms.
"""

from __future__ import division # 1/3 returns 0.333... instead of 0
from collections import defaultdict
from operator import add
from functools import reduce # Python 3 compability

from sympy import *
from sympy.physics import units # Sympy does not support
                                # PyPI package "Quantities"
from prj_helpers import adv_memoize
from water_permittivity import eps
from IAPWS95_density import get_water_density_derivatives, P_, T_
from manuscript_constants import P_, T_, P0, Tdash, eps0, alpha_LS, \
     M_W, N_A, gamma_prime_val, chi_prime, chi_tilde_minus_prime, \
     q_I_val, R_I_val


# Global variables (only LS scheme implemented)
COR_TYPES = ('B','C1','C2','D')
Y_TYPES   = ('G', 'H', 'S', 'CP', 'V', 'KT', 'AP')
Y_UNITS   = {'G':  units.joule/units.mole,
	     'H':  units.joule/units.mole,
	     'S':  units.joule/units.mole/units.kelvin,
	     'CP': units.joule/units.mole/units.kelvin,
	     'V':  units.joule/units.mole/units.pascal,
	     'KT': units.joule/units.mole/units.pascal**2,
	     'AP': units.joule/units.mole/units.pascal/units.kelvin
	     }

Y_UNITS_STR = {'G':  'J/mol',
	       'H':  'J/mol',
	       'S':  'J/(K mol)',
	       'CP': 'J/(K mol)',
	       'V':  'm^3/mol',
	       'KT': 'm^3/(Pa mol)',
	       'AP': 'm^3/(K mol)'
	     }
# Ions (only sodium and chloride ions implemented)
IONS = ('sod','cls')
ION_NAMES = {'sod': 'Na+', 'cls': 'Cl-'}

# Epsilon' is scaled down version of real water:
eps_prime = eps * 66.6/eps.subs({P_:P0, T_:Tdash})

rho_prime = symbols('rho_prime', cls=Function)(P_,T_)

q_I, R_I, N_W, gamma_prime = symbols('q_I, R_I, N_W, gamma_prime')


# Dependent parameters
# Computational box length as CALCULATED in Eq.37
# from density and number of water molecules
L=(N_W*M_W/N_A/rho_prime+4*pi/3*R_I**3)**sympify(1/3)

# Correction terms appropriate for Lattice summation
# Eq. 36 p. 17 in submitted manuscript
Delta_G_LS = {'B': 1/(8*pi*eps0)*N_A*q_I**2* \
	      (1-1/eps_prime)/L*(alpha_LS+4*pi/3*(R_I/L)**2 - \
				 16*pi**2/45*(R_I/L)**5),
	      'C1': -N_A/(6*eps0)*N_W*gamma_prime*q_I/L**3,
	      'C2': -N_A*q_I*4*pi*R_I**3/3/L**3 * \
			(chi_prime+chi_tilde_minus_prime/R_I),
	      'D':  1/(8*pi*eps0)*N_A*q_I**2 * \
				(1/eps-1/eps_prime)/R_I
	      }


# Delta_Y_LS is a dict of dicts, first key is `Y`, second key is cor_type:
Delta_Y_LS = defaultdict(dict)

for cor_type, cor_eq in Delta_G_LS.iteritems():
    Delta_Y_LS['G'][cor_type]  = cor_eq
    Delta_Y_LS['H'][cor_type]  = cor_eq - T_*diff(cor_eq,T_)
    Delta_Y_LS['S'][cor_type]  =            -diff(cor_eq,T_)
    Delta_Y_LS['CP'][cor_type] =         -T_*diff(cor_eq,T_,2)
    Delta_Y_LS['V'][cor_type]  =             diff(cor_eq,P_)
    Delta_Y_LS['KT'][cor_type] =            -diff(cor_eq,P_,2)
    Delta_Y_LS['AP'][cor_type] =             diff(cor_eq,P_,T_)

def get_rho_subs(P_val, T_val, reltol=None, verbose=False,
		 IAPWS95_verbose=False):
    from functools import reduce
    from operator import add

    """
    The IAPWS95 expression is implicit and rho must
    be calculated numerically
    """
    if reltol == None: reltol = 1e-9
    scaling_factor = 968.2/997.047625482 # See Table 3 on p. 43 in MS
    if verbose: print "Calculating water density P and T derivatives using IAPWS95..."
    P_order = 2 # We need up to second derivative wrt to P
    T_order = 2 # We need up to second derivative wrt to T
    skip_sigs = (((P_, 1), (T_, 2)),
		 ((P_, 2), (T_, 2)),
		 ((P_, 2), (T_, 1)))

    val, err = get_water_density_derivatives(P_order, T_order, P_val,
					     T_val, None,
					     IAPWS95_verbose,
					     reltol,
					     skip_sigs=skip_sigs)

    subs = {}
    for k, v in val.iteritems():
	subs[Derivative(rho_prime, *reduce(add, k))] =v*scaling_factor
	if verbose:
	    k = Derivative(rho_prime, *reduce(add, k))
	    print k,'=',subs[k]
    return subs


def test_get_rho_subs():
    print get_rho_subs(P0, Tdash)


@adv_memoize()
def get_Delta_Y_cor_LS(Y, P_val, T_val, N_W_val, ion, cor_type="all", verbose=False):
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
    assert(ion in IONS)
    assert(Y in Y_TYPES)

    try:
	float(P_val/units.pascal)
    except:
	# Assume Pascal
	P_val *= units.pascal

    try:
	float(T_val/units.kelvin)
    except:
	# Assume Kelvin
	T_val *= units.kelvin

    subsd = get_rho_subs(P_val, T_val, None, verbose)
    if cor_type == "all":
	if verbose: print "Calculating all correction term types ("+\
	   ", ".join(COR_TYPES)+")"
	result = {}
	for key in Delta_Y_LS[Y].keys():
	    # Call the the function itself for each Y
	    result[key] = get_Delta_Y_cor_LS(Y, P_val, T_val,
					     N_W_val, ion, key,
					     verbose)[key]
	return result
    else:
	if verbose: print "Calculating cor_type: {}".format(cor_type)
	subsd.update({P_: P_val, T_: T_val, N_W: N_W_val,
		 q_I: q_I_val[ion], R_I: R_I_val[ion],
		 gamma_prime: gamma_prime_val})
	return {cor_type: Delta_Y_LS[Y][cor_type].subs(subsd).evalf().simplify()}


def test_get_Delta_Y_cor_LS(verbose=False):
    Y='G'; P_val=P0; T_val = Tdash; N_W_val=1024; ion='sod'
    print get_Delta_Y_cor_LS(Y, P_val,T_val,N_W_val,ion,"all",
			     verbose)

