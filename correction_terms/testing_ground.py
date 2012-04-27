#!/usr/bin/env python

from IAPWS_density import *

# print a.evalf()
# solve(pressure_relation, rho)


# get_water_density(P0,Tdash)

rho0 = sympify(1000.0)*units.kg/units.meter**3
pressure_relation.subs({P:P0,T:Tdash,rho:rho0}).evalf()
