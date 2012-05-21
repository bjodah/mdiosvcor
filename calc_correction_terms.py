#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This work is open source and is released under the
# 2-clause BSD license (see LICENSE.txt for further information)
# Copyright (c) 2011, 2012, Björn Ingvar Dahlgren

"""
Program for calculating the correction terms applied in
the publication:

   Björn Dahlgren, Maria M. Reif, Philippe H. Hünenberger & Niels Hansen
   "Calculation of derivative thermodynamic hydration and aqueous partial
   molar properties of ions based on atomistic simulations"
   J. Chem. Theory Comput.

The correction terms are to be applied to thermodynamic solvation parameters and partial molar variables determined by MD simulation.

Latest version is available at github.com/bjodah/mdiosvcor

Implemented for use in research project in the IGC group at ETH

If you make use the code in your own research, please
cite the article referenced above.
See README.md and LICENSE.txt for further information
"""
import argparse, os, sys

absdirname = os.path.abspath(os.path.dirname(sys.argv[0]))
os.environ['MEMOIZE_CACHE_DIR'] = os.path.join(absdirname,'cache/')

from mdiosvcor.manuscript_equations import get_Delta_Y_cor_LS, COR_TYPES, Y_TYPES, IONS, ION_NAMES

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

    print(get_Delta_Y_cor_LS(Y	      = argd['property'],
			     P_val    = argd['pressure'],
			     T_val    = argd['temperature'],
			     N_W_val  = argd['nwater'],
			     ion      = argd['ion'],
			     cor_type = argd['cortype'],
			     verbose  = argd['verbose']))
