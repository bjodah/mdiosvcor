#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This work is open source and is released under the
# 2-clause BSD license (see LICENSE.txt for further information)
# Copyright (c) 2011, 2012, Björn Ingvar Dahlgren

import sympy.core
from sympy import *
from sympy.physics import units

from manuscript_constants import P0, Tdash

"""
Permitivity of water according to
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

"""
# Variables
P_  = Symbol('P') # Pressure (Intensive state variable)
T_  = Symbol('T') # Temperature (Intensive state variable)

U = [0,
  3.4279e2,
 -5.0866e-3 / units.kelvin,
  9.4690e-7 / units.kelvin / units.kelvin,
 -2.0525,
  3.1159e3  * units.kelvin,
 -1.8289e2  * units.kelvin,
 -8.0325e3  * units.bar,
  4.2142e6  * units.kelvin * units.bar,
  2.1417    / units.kelvin * units.bar]
B  = U[7] + U[8]/T_ + U[9]*T_
C  = U[4] + U[5]/(U[6]+T_)
eps1000 = U[1]*exp(U[2]*T_+U[3]*T_**2)
eps = eps1000 + C*ln((B+P_)/(B+1000.0*units.bar))

def get_water_eps(P_val,T_val, expr=eps):
    # If P is without unit, assume Pascal:
    try:
        float(P_val/units.pascal)
    except:
        P_val *= units.pascal

    # If T is without unit, assume Kelvin:
    try:
        float(T_val/units.kelvin)
    except:
        T_val *= units.kelvin

    return float(expr.subs({P_:P_val, T_:T_val}))
