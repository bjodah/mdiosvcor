# -*- coding: utf-8 -*-
# <nbformat>2</nbformat>

# <codecell>

#!/usr/bin/env python
import sympy.core
from sympy import *
from sympy.physics import units

# <codecell>

# Permitivity of water according Bradley and Pitzer

# Variables
P  = Symbol('P') # Pressure (Intensive state variable)
T  = Symbol('T') # Temperature (Intensive state variable)

U1 =  3.4279e2
U2 = -5.0866e-3 / units.kelvin
U3 =  9.4690e-7 / units.kelvin / units.kelvin
U4 = -2.0525
U5 =  3.1159e3  * units.kelvin
U6 = -1.8289e2  * units.kelvin
U7 = -8.0325e3  * units.bar
U8 =  4.2142e6  * units.kelvin * units.bar
U9 =  2.1417    / units.kelvin * units.bar
B  = U7 + U8/T + U9*T
C  = U4 + U5/(U6+T)
eps1000 = U1*exp(U2*T+U3*T**2)
eps = eps1000 + C*ln((B+P)/(B+1000.0*units.bar))

eps_prime = eps * 66.6/78.4

# <codecell>

def get_eps(P_val,T_val):
    # If P is without unit:
    try:
        float(P_val/units.pascal)
    except:
        unit = units.pascal
        print 'Warning: assumed P_val was in unit:' + str(unit)
        P_val *= unit
    
    # If T is without unit
    try:
        float(T_val/units.kelvin)
    except:
        unit = units.kelvin
        print 'Warning: assumed T_val was in unit:' + str(unit)
        T_val *= unit
        
    return float(eps.subs({P:P_val,T:T_val}))

def get_eps_prime(P_val,T_val):
    # If P is without unit:
    try:
        float(P_val/units.pascal)
    except:
        unit = units.pascal
        print 'Warning: assumed P_val was in unit:' + str(unit)
        P_val *= unit
    
    # If T is without unit
    try:
        float(T_val/units.kelvin)
    except:
        unit = units.kelvin
        print 'Warning: assumed T_val was in unit:' + str(unit)
        T_val *= unit
    
    return float(eps_prime.subs({'P':P_val,'T':T_val}))

# <codecell>

get_eps(101.3e3,298.15)

