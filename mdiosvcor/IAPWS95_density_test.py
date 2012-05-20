#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This work is open source and is released under the
# 2-clause BSD license (see LICENSE.txt for further information)
# Copyright (c) 2011, 2012, Björn Ingvar Dahlgren

"""
Unit-tests for IAPWS95_density.py
"""

import unittest
import numpy as np
#import numpy as np
from sympy import *
from IAPWS95_density import *
from helpers import get_sympified


# From Table 6.6 p. 436 (50 in PDF)

phi0_ref_500 = {
                 'f'             : 0.204797733e1,
                 'dfddelta'      : 0.384236747,
                 'd2fddelta2'    :-0.147637878,
                 'dfdtau'        : 0.904611106e1,
                 'd2fdtau2'      :-0.193249185e1,
                 'd2fddeltadtau' : 0,
                }

phi__r_ref_500 = {
                 'f'             :-0.342693206e1,
                 'dfddelta'      :-0.364366650,
                 'd2fddelta2'    : 0.856063701,
                 'dfdtau'        :-0.581403435e1,
                 'd2fdtau2'      :-0.223440737e1,
                 'd2fddeltadtau' :-0.112176915e1,
                 }

state_500_subs = {tau_: val_T_c / (sympify(500.0) * units.kelvin),
                  delta_: sympify(838.025) * units.kg / units.meter**3\
		         / val_rho_c
                  }

phi0_ref_647 =  {
                 'f'             :-0.156319605e1,
                 'dfddelta'      : 0.899441341,
                 'd2fddelta2'    :-0.808994726,
                 'dfdtau'        : 0.980343918e1,
                 'd2fdtau2'      :-0.343316334e1,
                 'd2fddeltadtau' : 0.0,
                }

phi__r_ref_647 = {
                  'f'             :-0.121202657e1,
                  'dfddelta'      :-0.714012024,
                  'd2fddelta2'    : 0.475730696,
                  'dfdtau'        :-0.321722501e1,
                  'd2fdtau2'      :-0.996029507e1,
                  'd2fddeltadtau' :-0.133214720e1,
                 }

state_647_subs = {tau_:   val_T_c / (sympify(647.0) * units.kelvin),
                  delta_: sympify(358.0) * units.kg / units.meter**3 \
		         / val_rho_c
                  }

global derivative_order
derivative_order = ('f','dfddelta','d2fddelta2','dfdtau','d2fdtau2','d2fddeltadtau')


derivative_vars = {
    # derivative     (der_tuple, factor)
    # der_tuple: what vairables to derivate with respect to
    # factor what to multiply resulting derivative with in
    # order to get the expected derivative.
    # Note: The reason why the derivation is not done
    #       explicitly with respect to tau/delta respectively
    #       is that sympy complains on rho_c/T_c having units
    'f'             :None,
    'dfddelta'      :(delta_, 1),
    'd2fddelta2'    :(delta_, 2),
    'dfdtau'        :(tau_, 1),
    'd2fdtau2'      :(tau_, 2),
    'd2fddeltadtau' :(delta_, 1, tau_, 1)
}

refs = [(state_500_subs, phi0_,   phi0_ref_500),
        (state_500_subs, phi__r_, phi__r_ref_500),
	(state_647_subs, phi0_,   phi0_ref_647),
	(state_647_subs, phi__r_, phi__r_ref_647)]

class Test_ref(unittest.TestCase):
    """
    Verify that the verfication data is reproduced
    """
    derivative_order = derivative_order
    derivative_vars  = derivative_vars

    def test_verification_data(self):
	for subs_d,func,ref in refs:
	    ref_vals = np.array([ref[key] for key in self.derivative_order])
	    calc_vals = []
	    for var in self.derivative_order:
		wrt = self.derivative_vars[var]
		if wrt != None:
		    deriv = func.diff(*wrt)
		else:
		    deriv = func
		calc_vals.append(deriv.subs(subs_d))
	    calc_vals = np.array([float(x) for x in calc_vals])
	    self.assertTrue(np.allclose(calc_vals,ref_vals),
			    'Discrepancy in %s' % str(subs_d))


if __name__ == '__main__':
    unittest.main()
