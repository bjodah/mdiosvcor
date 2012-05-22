#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This work is open source and is released under the
# 2-clause BSD license (see LICENSE.txt for further information)
# Copyright (c) 2011, 2012, Björn Ingvar Dahlgren

"""
Program for calculating the correction terms applied in
the publication:

   Björn Dahlgren, Maria M. Reif, Philippe H. Hünenberger & Niels Hansen
   `Calculation of derivative thermodynamic hydration and aqueous partial
   molar properties of ions based on atomistic simulations`
   J. Chem. Theory Comput.

The correction terms are to be applied to thermodynamic solvation parameters and partial molar variables determined by MD simulation.

Latest version is available at github.com/bjodah/mdiosvcor

Implemented for use in research project in the IGC group at ETH

If you make use the code in your own research, please
cite the article referenced above.
See README.md and LICENSE.txt for further information
"""
import argparse, os, sys

if __name__ == '__main__':
    # Running as script
    absdirname = os.path.abspath(os.path.dirname(sys.argv[0]))
else:
    # Being imported
    absdirname = os.getcwd()

os.environ['MEMOIZE_CACHE_DIR'] = os.path.join(absdirname,'cache/')

from mdiosvcor.manuscript_equations import get_Delta_Y_cor_LS, COR_TYPES, Y_TYPES, Y_UNITS, Y_UNITS_STR, IONS, ION_NAMES
from mdiosvcor.prj_helpers import get_unitless

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, epilog="Abbrevations: "+\
				     "G = Gibb's free energy, H = Enthalpy, S = Entropy, "+\
				     "CP = Isobaric heat capacity, V = molar volume, "+\
				     "KT = isothermal volume-compressibility, AP = isobaric volume-expansivity")
    parser.add_argument('-P', '--pressure',	type=float,	default=101300, help='Pressure in Pascal')
    parser.add_argument('-T', '--temperature',	type=float,     default=298.15, help='Temperature in Kelvin')
    parser.add_argument('-N', '--nwater',	type=int,       default=1024,   help='Number of water molecules')
    parser.add_argument('-Y', '--property',	type=str,	default='G',    help='Thermodynamic property, specify one of: '+', '.join(Y_TYPES))
    parser.add_argument('-c', '--cortype',	type=str,	default='all',  help='Correction type, specify either \'all\' or any of: '+', '.join(COR_TYPES))
    parser.add_argument('-i', '--ion',		type=str,	default='sod',  help='Ion type, keywords: '+', '.join([k+": "+v for k,v in ION_NAMES.items()]))
    parser.add_argument('-v', '--verbose',    action='store_true',              help='Enable verbose output.')

    args = parser.parse_args()
    argd = vars(args)
    if argd['verbose']:
	fstr = "Calculating for P={}, T={}, N_W_val={}, ion={}"
	print fstr.format(*[argd[k] for k in ('pressure',
					      'temperature',
					      'nwater',
					      'ion')])

    result = get_Delta_Y_cor_LS(argd['property'],
				argd['pressure'],
				argd['temperature'],
				argd['nwater'],
				argd['ion'],
				argd['cortype'],
				argd['verbose'])

    Y = argd['property']
    try:
	for key in COR_TYPES:
	    fmtstr = "{0: >2}: {1} {2}"
	    # get_unitless is invoked due to residual m**2.0/m**2
	    print fmtstr.format(key, str(get_unitless(result[key]/ \
						      Y_UNITS[Y])),
				Y_UNITS_STR[Y])
    except:
	raise


def batch_calc(Ys, Ps, Ts, NWs, Is, cors, verbose=False, dump_to_file=None):
    from itertools import product
    result = {}
    fmtstr="Y={}, P={}, T={}, NW={}, I={}, COR={}: {} {}"
    for conditions in product(Ys, Ps, Ts, NWs, Is, cors):
	result[conditions] = get_Delta_Y_cor_LS(*(conditions+(verbose,)))
	Y   = conditions[0]
	cor = conditions[5]
	unitless_val = get_unitless(result[conditions][cor]/Y_UNITS[Y])
	str_unit = Y_UNITS_STR[Y]
	print fmtstr.format(*(conditions+(unitless_val, str_unit)))

    if dump_to_file:
	try:
	    import cPickle as pickle
	except ImportError:
	    import pickle
	pickle.dump(result, open(dump_to_file, 'wb'))


def test_batch_calc():
    batch_calc(['G','H'],[101.3e3,101.3e3*5e3],[298.15],[1024],['sod'],['B','C1',],verbose=True,dump_to_file='ignore_pickle_dump_test_batch_calc')


def calc_all_for_project():
    batch_calc(Y_TYPES,[101.3e3],
	       [273.15, 285.65, 298.15, 310.65, 323.15],
	       [1024, 724, 512], ['sod','cls'],COR_TYPES,
	       dump_to_file="ignore_all_cor_prj_T")
    batch_calc(Y_TYPES,[1e5, 1e5*5e3, 1e5*10e3], [298.15],
	       [1024, 724, 512], ['sod','cls'],COR_TYPES,
	       dump_to_file="ignore_all_cor_prj_P")
