#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This work is open source and is released under the
# 2-clause BSD license (see LICENSE.txt for further information)
# Copyright (c) 2011, 2012, Björn Ingvar Dahlgren

"""
Program for calculating permitivity of water according to
the parametrization by Bradley and Pitzer

         Bradley, D.J.; Pitzer, K.S. `Thermodynamics of electrolytes. 12. Dielectric properties of water and Debye--Hueckel parameters to 350/sup 0/C and 1 kbar`, J. Phys. Chem.; Journal Volume 83 (12) (1979), pp. 1599-1603
         http://pubs.acs.org/doi/abs/10.1021/j100475a009
         DOI: 10.1021/j100475a009

See README.md and LICENSE.txt for further information
"""
import argparse, os, sys

from mdiosvcor.water_permittivity import eps, P_, T_
from sympy import Derivative, symbols, Function
from sympy.physics import units
from mdiosvcor.prj_helpers import get_unitless

absdirname = os.path.abspath(os.path.dirname(sys.argv[0]))
if not 'MEMOIZE_CACHE_DIR' in os.environ:
    os.environ['MEMOIZE_CACHE_DIR'] = os.path.join(absdirname,'cache/')


def get_water_eps_derivatives(P_order, T_order,
			      val_P=None,
			      val_T=None,
			      expr=eps,
			      ret_w_units=True
			      ):
    """
    Function to calculate arbitrary P/T derivative of water's
    relative permittivity according to parametrization by
    Bradley and Pitzer

    :Optional arguments:
    expr
        Expression of eps as function of P, T (sympy symbols)
	Default is the parametrization by Bradley and Pitzer
    ret_w_units
        Should the calculated values be returned with or without
	units (sympy.phiscs.units). Default is True
    """
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

    expr_deriv = expr.diff(P_,P_order,T_,T_order)
    result = expr_deriv.subs({P_: val_P, T_: val_T})
    if not ret_w_units:
	return get_unitless(result)
    else:
	return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-T', '--temperature', type=float,  default=298.15, help='Temperature in Kelvin')
    parser.add_argument('-P', '--pressure',    type=float,  default=101300, help='Pressure in Pascal')
    parser.add_argument('-i', '--porder',      type=int,    default=0,      help='Derivative order with respect to presusre')
    parser.add_argument('-j', '--torder',      type=int,    default=0,      help='Derivative order with respect to temperature')
    parser.add_argument('-u', '--units',     action='store_true',           help='Print out units of calculated density (derivative).')
    parser.add_argument('-v', '--verbose',   action='store_true',           help='Print numerical convergence info etc.')

    args = parser.parse_args()
    argd = vars(args)

    calc_prop = Derivative(symbols('eps',cls=Function),
			   P_, argd['porder'], T_, argd['torder'])

    if argd['verbose']:
	fstr = "Calculating {} for water (eps = relative permitivity) at P = {} Pa, T = {} K"
	print fstr.format(calc_prop, *[argd[k] for k in ('pressure', 'temperature')]
			  )
    val = get_water_eps_derivatives(argd['porder'],
				    argd['torder'],
				    argd['pressure'],
				    argd['temperature'],
				    eps,
				    argd['units'])
    print calc_prop,':', val
