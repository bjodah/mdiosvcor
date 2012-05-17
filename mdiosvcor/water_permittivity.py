#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sympy.core
from sympy import *
from sympy.physics import units

from manuscript_constants import P0, Tdash

# Permitivity of water according Bradley and Pitzer

# Variables
P  = Symbol('P') # Pressure (Intensive state variable)
T  = Symbol('T') # Temperature (Intensive state variable)

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
B  = U[7] + U[8]/T + U[9]*T
C  = U[4] + U[5]/(U[6]+T)
eps1000 = U[1]*exp(U[2]*T+U[3]*T**2)
eps = eps1000 + C*ln((B+P)/(B+1000.0*units.bar))

def get_water_eps(P_val,T_val, expr=eps):
    # If P is without unit, assume Pascal:
    try:
        float(P_val/units.pascal)
    except:
        unit = units.pascal
        print 'Warning: assumed P_val was in unit:' + str(unit)
        P_val *= unit

    # If T is without unit, assume Kelvin:
    try:
        float(T_val/units.kelvin)
    except:
        unit = units.kelvin
        print 'Warning: assumed T_val was in unit:' + str(unit)
        T_val *= unit

    return float(expr.subs({P:P_val,T:T_val}))
