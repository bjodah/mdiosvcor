#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This work is open source and is released under the
# 2-clause BSD license (see LICENSE.txt for further information)

from __future__ import division     # 1/3 returns 0.333... instead of 0
from collections import defaultdict
from functools import reduce        # For Python 3 compability
from operator import add

from sympy import *
from sympy.physics import units # Sympy does not support
                                # PyPI package "Quantities"

from water_permittivity import eps, get_water_eps
from IAPWS95_density import get_water_density_derivatives
from manuscript_constants import P, T, P0, Tdash, eps0, alpha_LS, \
     M_W, N_A, gamma_prime_val, chi_prime, chi_tilde_minus_prime, \
     q_I_val, R_I_val


# Global variables (only LS scheme implemented)
COR_TYPES = ('B','C1','C2','D')
Y_TYPES   = ('G', 'H', 'S', 'CP', 'V', 'KT', 'AP')
# Ions (only sodium and chloride ions implemented)
IONS = ('sod','cls')
ION_NAMES = {'sod': 'Na+', 'cls': 'Cl-'}

# Epsilon' is scaled down version of real water:
eps_prime = eps * 66.6/get_water_eps(P0,Tdash)

rho_prime = symbols('rho_prime', cls=Function)(P,T)

q_I, R_I, N_W, gamma_prime = symbols('q_I, R_I, N_W, gamma_prime')


# Dependent parameters
# Computational box length as CALCULATED in Eq.37
# from density and number of water molecules
L=(N_W*M_W/N_A/rho_prime+4*pi/3*R_I**3)**(1/3)

# Correction terms appropriate for Lattice summation
# Eq. 36 p. 17 in submitted manuscript
Delta_G_LS = {'B':  (8*pi*eps0)*N_A*q_I**2* \
	      (1-1/eps_prime)/L*(alpha_LS+4*pi/3*(R_I/L)**2 - \
				 16*pi**2/45*(R_I/L)**5),
	      'C1': -N_A/(6*eps0)*N_W*gamma_prime*q_I/L**3,
	      'C2': -N_A*q_I*4*pi*R_I**3/3/L**3 * \
			(chi_prime+chi_tilde_minus_prime/R_I),
	      'D':  1/8/pi/eps0*N_A*q_I**2 * \
				(1/eps-1/eps_prime)/R_I
	      }


# Delta_Y_LS is a dict of dicts, first key is `Y`, second key is cor_type:
Delta_Y_LS	= defaultdict(dict)

for cor_type, cor_eq in Delta_G_LS.iteritems():
    Delta_Y_LS['G'][cor_type]  = cor_eq
    Delta_Y_LS['H'][cor_type]  = 1 - T*diff(cor_eq,T)
    Delta_Y_LS['S'][cor_type]  = -diff(cor_eq,T)
    Delta_Y_LS['CP'][cor_type] = -T*diff(cor_eq,T,2)
    Delta_Y_LS['V'][cor_type]  = diff(cor_eq,P)
    Delta_Y_LS['KT'][cor_type] = -diff(cor_eq,P,2)
    Delta_Y_LS['AP'][cor_type] = diff(cor_eq,P,T)


def get_rho_subs(P_val, T_val, abstol=1e-9, verbose=False):
    """
    The IAPWS95 expression is implicit and rho must
    be calculated numerically
    """
    scaling_factor = 968.2/997.05 # See Table 3 on p. 43 in MS
    if verbose: print "Calculating water density P and T derivatives using IAPWS95..."
    P_order = 2 # We need up to second derivative wrt to P
    T_order = 2 # We need up to second derivative wrt to T
    val, err = get_water_density_derivatives(P_order, T_order, P_val,
					     T_val, None, verbose, abstol)
    subs = {}
    for k,v in val:
	subs[Derivative(rho_prime, *reduce(add, k))] = v * scaling_factor
    return subs


def test_get_rho_subs():
    print get_rho_subs(P0, Tdash)


def get_Delta_Y_cor_LS(Y, P_val, T_val, N_W_val, ion,cor_type="all",verbose=False):
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

    if verbose: print "Calculating for P={}, T={}".format(P_val,T_val)
    if cor_type == "all":
	if verbose: print "Calculating all correction term types ("+", ".join(COR_TYPES)+")"
	result = {}
	for key in Delta_Y_LS[Y].keys():
	    result[key] = get_Delta_Y_cor_LS(Y, P_val, T_val,N_W_val,
					     ion, cor_type=key, verbose=verbose)
	return result
    else:
	if verbose: print "Calculating cor_type: {}".format(cor_type)
	subsd = {P: P_val, T: T_val, N_W: N_W_val,
		 q_I: q_I_val[ion], R_I: R_I_val[ion],
		 gamma_prime: gamma_prime_val}
	subsd.update(get_rho_subs(P_val, T_val, verbose=verbose))
	return Delta_Y_LS[Y][cor_type].subs(subsd)


def test_get_Delta_Y_cor_LS(verbose=False):
    Y='G'; P_val=P0; T_val = Tdash; N_W_val=1024; ion='sod'
    print get_Delta_Y_cor_LS(Y, P_val,T_val,N_W_val,ion,"all",
			     verbose=verbose)

