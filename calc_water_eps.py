#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This work is open source and is released under the
# 2-clause BSD license (see LICENSE.txt for further information)
# Copyright (c) 2011, 2012, Björn Ingvar Dahlgren

"""
Program for calculating relative permitivity of water according to
the parametrization by Bradley and Pitzer

         Bradley, D.J.
         Pitzer, K.S.
         Journal Name: J. Phys. Chem.; (United States); Journal Volume: 83:12
         1979
         J. Phys. Chem.; (United States); Journal Volume: 83:12
         Pages: 1599-1603
         Thermodynamics of electrolytes. 12. Dielectric properties of water and Debye--Hueckel parameters to 350/sup 0/C and 1 kbar
         http://pubs.acs.org/doi/abs/10.1021/j100475a009
         DOI: 10.1021/j100475a009

See README.md and LICENSE.txt for further information
"""
import argparse, os, sys

absdirname = os.path.abspath(os.path.dirname(sys.argv[0]))
os.environ['MEMOIZE_CACHE_DIR'] = os.path.join(absdirname,'cache/')

from mdiosvcor.water_permittivity import get_water_eps

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-T', '--temperature',	type=float,     default=298.15, help='Temperature in Kelvin')
    parser.add_argument('-P', '--pressure',	type=float,	default=101300, help='Pressure in Pascal')
    parser.add_argument('-v', '--verbose',    action='store_true',              help='Print numerical convergence info etc.')

    args = parser.parse_args()
    argd = vars(args)
    if argd['verbose']:
	fstr = "Calculating relative permitivity of water for P = {} Pa, T = {} K"
	print fstr.format(*[argd[k] for k in ('pressure',
					      'temperature')
					      ])
    val = get_water_eps( P_val    = argd['pressure'],
			 T_val    = argd['temperature'])
    print "val:   ", val
