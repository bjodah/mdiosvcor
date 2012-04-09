#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
import numpy as np
#import numpy as np
from sympy import *
from IAPWS95_density import *

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

state_500_subs = {T:   500*units.kelvin,
                  rho: 838.025 * units.kg / units.meter**3
                  }

derivative_order = ('f','dfddelta','d2fddelta2','dfdtau','d2fdtau2','d2fddeltadtau')

derivative_vars = {
    # derivative     (der_tuple, factor)
    # der_tuple: what vairables to derivate with respect to
    # factor what to multiply resulting derivative with in
    # order to get the expected derivative.
    # Note: The reason why the derivation is not done
    #       explicitly with respect to tau/delta respectively
    #       is that sympy complains on rho_c/T_c having units
    'f'             :((None,),  1),
    'dfddelta'      :((rho,),  rho_c),
    'd2fddelta2'    :((rho,2), rho_c**2),
    'dfdtau'        :((T,),    T_c),
    'd2fdtau2'      :((T,2),   T_c**2),
    'd2fddeltadtau' :((rho,T), rho_c*T_c),
}

refs = [(state_500_subs, phi0,   phi0_ref_500),
        (state_500_subs, phi__r, phi__r_ref_500)]

class Test_ref(unittest.TestCase):
    """
    Verify that the verfication data is reproduced
    """

    def runTest(self):
        test_verification_data()
        
    def test_verification_data(self):
        for subs,func,ref in refs:
            ref_vals = np.array([ref[key] for key in derivative_order])
            calc_vals = []
            for var in derivate_order:
                wrt, factor = derivative_vars[var]
                expr = func.diff(func,*wrt)*factor
                calc_vals.append(expr.subs(subs))
            self.assertTrue(np.allclose(calc_vals,ref_vals))


