#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This work is open source and is released under the
# 2-clause BSD license (see LICENSE.txt for further information)
# Copyright (c) 2011, 2012, Björn Ingvar Dahlgren

"""
Program for calculating the density of water
as function pressure and temperature
according to IAPWS95 formulation.

The code is entirely based on the publication:
   Wagner, W., and A. Pruß. “The IAPWS Formulation 1995 for the
   Thermodynamic Properties of Ordinary Water Substance for
   General and Scientific Use.” Journal of Physical and Chemical
   Reference Data 31, no. 2 (June 7, 2002): 387–535.

Latest version is available at github.com/bjodah/mdiosvcor

Implemented for use in a research project in the IGC group at ETH.

See README.md and LICENSE.txt for further information
"""
import argparse, os, sys

absdirname = os.path.abspath(os.path.dirname(sys.argv[0]))
os.environ['MEMOIZE_CACHE_DIR'] = os.path.join(absdirname,'cache/')

from mdiosvcor.IAPWS95_density import (
    get_water_density_derivatives, P_, T_
    )

from sympy import Derivative, symbols, Function

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-P', '--pressure',    type=float, default=101300, help='Pressure in Pascal')
    parser.add_argument('-T', '--temperature', type=float, default=298.15, help='Temperature in Kelvin')
    parser.add_argument('-i', '--porder',      type=int,   default=0,      help='Derivative order with respect to presusre')
    parser.add_argument('-j', '--torder',      type=int,   default=0,      help='Derivative order with respect to temperature')
    parser.add_argument('-r', '--reltol',      type=float, default=1e-9,   help='Relative tolerance for root finding algorithm.')
    parser.add_argument('-u', '--units',     action='store_true',          help='Print out units of calculated density (derivative).')
    parser.add_argument('-v', '--verbose',   action='store_true',          help='Print numerical convergence info etc.')

    args = parser.parse_args()
    argd = vars(args)

    calc_prop = Derivative(symbols('rho',cls=Function),
			   P_, argd['porder'], T_, argd['torder'])

    if argd['verbose']:
	fstr = "Calculating {} (rho = water density) at P={}, T={}, with a relative tolerance of {}"
	print fstr.format(calc_prop, *[argd[k] for k in ('pressure',
					      'temperature',
					      'reltol')])

    use_numexpr = False; use_finite_difference = False
    val, err = get_water_density_derivatives(argd['porder'],
					argd['torder'],
					argd['pressure'],
					argd['temperature'],
					None,
					argd['verbose'],
					argd['reltol'],
					use_numexpr,
					use_finite_difference,
					argd['units'])

    print "{}: {}".format(calc_prop,
	val[((P_,argd['porder']),(T_,argd['torder']))])
    print "trunc_err: {}".format(
	abs(err[((P_,argd['porder']),(T_,argd['torder']))]))
