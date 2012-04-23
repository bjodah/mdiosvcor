#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Provides constants used in the IAPWS95 formulation.

# Author: Björn Dahlgren
# Implemented for use in research project in the IGC group at ETH

# To the extent possible under law, Bjoern Dahlgren has waived all
# copyright and related or neighboring rights to this work.

from __future__ import division
from sympy import *
from sympy.physics import units
from root_finding import find_root

# Variables
P  = Symbol('P')   # Pressure (Intensive state variable)
T  = Symbol('T')   # Temperature (Intensive state variable)
rho= Symbol('rho') # Denisty (Intensive state variable)


# Constants used in the master equation

_n0 = [sympify(x) for x in \
        [0.0,
       -8.3204464837497, 6.6832105275932, 3.00632,
        0.012436,        0.97315,         1.27950,
        0.96956,         0.24873]]


_gamma0 = [sympify(x) for x in \
           [0.0,
           0.0,               0.0,               0.0,
           1.28728967,      3.53734222,      7.74073708,
           9.24437796,     27.5075105]]


# The raw data is copied from IAPWS95.py by Kiran Pashikanti
_n =  [sympify(x) for x in \
       [0.0,
          0.12533547935523e-1,  0.78957634722828e+1, -0.87803203303561e+1,
          0.31802509345418e+0, -0.26145533859358e+0, -0.78199751687981e-2,
          0.88089493102134e-2, -0.66856572307965e+0,  0.20433810950965e+0,
         -0.66212605039687e-4, -0.19232721156002e+0, -0.25709043003438e+0,
          0.16074868486251e+0, -0.40092828925807e-1,  0.39343422603254e-6,
         -0.75941377088144e-5,  0.56250979351888e-3, -0.15608652257135e-4,
          0.11537996422951e-8,  0.36582165144204e-6, -0.13251180074668e-11,
         -0.62639586912454e-9, -0.10793600908932e+0,  0.17611491008752e-1,
          0.22132295167546e+0, -0.40247669763528e+0,  0.58083399985759e+0,
          0.49969146990806e-2, -0.31358700712549e-1, -0.74315929710341e+0,
          0.47807329915480e+0,  0.20527940895948e-1, -0.13636435110343e+0,
          0.14180634400617e-1,  0.83326504880713e-2, -0.29052336009585e-1,
          0.38615085574206e-1, -0.20393486513704e-1, -0.16554050063734e-2,
          0.19955571979541e-2,  0.15870308324157e-3, -0.16388568342530e-4,
          0.43613615723811e-1,  0.34994005463765e-1, -0.76788197844621e-1,
          0.22446277332006e-1, -0.62689710414685e-4, -0.55711118565645e-9,
         -0.19905718354408e+0,  0.31777497330738e+0, -0.11841182425981e+0,
         -0.31306260323435e+2,  0.31546140237781e+2, -0.25213154341695e+4,
         -0.14874640856724e+0,  0.31806110878444e+0]]

_c =  [sympify(x) for x in \
       [0,
          0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
          2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
          2, 2, 2, 3, 3, 3, 3, 4, 6, 6, 6, 6]]
_d =  [sympify(x) for x in \
       [0,
          1, 1, 1, 2,  2,  3,  4,  1,  1, 1, 2,  2,  3,  4,
          4, 5, 7, 9, 10, 11, 13, 15,  1, 2, 2,  2,  3,  4,
          4, 4, 5, 6,  6,  7,  9,  9,  9, 9, 9, 10, 10, 12,
          3, 4, 4, 5, 14,  3,  6,  6,  6, 3, 3,  3]]
_t =  [sympify(x) for x in \
       [0,
          -0.5, 0.875,  1.0,  0.5,  0.75, 0.375,  1.0,  4.0,  6.0, 12.0,  1.0,
           5.0, 4.0  ,  2.0, 13.0,  9.0 , 3.0  ,  4.0, 11.0,  4.0, 13.0,  1.0,
           7.0, 1.0  ,  9.0, 10.0, 10.0 , 3.0  ,  7.0, 10.0, 10.0,  6.0, 10.0,
          10.0, 1.0  ,  2.0,  3.0,  4.0 , 8.0  ,  6.0,  9.0,  8.0, 16.0, 22.0,
          23.0,23.0  , 10.0, 50.0, 44.0, 46.0  , 50.0,  0.0,  1.0,  4.0]]

#offset constants, no point in creating giant nearly empty lists
_alpha = dict([(k,sympify(v)) for k,v in \
               {52:20.0,53:20.0,54:20.0}.iteritems()])
_beta  = dict([(k,sympify(v)) for k,v in \
               {52:150.0,53:150.0,54:250.0,55:0.3,56:0.3}.iteritems()])

_gamma = dict([(k,sympify(v)) for k,v in \
               {52:1.21,53:1.21,54:1.25}.iteritems()])
_epsilon = dict([(k,sympify(v)) for k,v in \
               {52:1.0, 53:1.0 ,54:1.0 }.iteritems()])
_a     = dict([(k,sympify(v)) for k,v in \
               {55:3.5 ,56:3.5}.iteritems()])
_b     = dict([(k,sympify(v)) for k,v in \
               {55:0.85,56:0.95}.iteritems()])
_B     = dict([(k,sympify(v)) for k,v in \
               {55:0.2,56:0.2}.iteritems()])
_C     = dict([(k,sympify(v)) for k,v in \
               {55:28.0,56:32.0}.iteritems()])
_D     = dict([(k,sympify(v)) for k,v in \
               {55:700.0,56:800.0}.iteritems()])
_A     = dict([(k,sympify(v)) for k,v in \
               {55:0.32,56:0.32}.iteritems()])

# <codecell>

R = sympify(461.51805) * units.joule / units.kelvin / units.kg

rho_c = sympify(322.0) * units.kg / units.meter**3
delta = rho/rho_c

T_c = sympify(647.096) * units.kelvin
tau = T_c/T

# P_c = 22.064e6 * units.pascal # Critical pressure, not used explicitly

# phi = f/RT = phi0 + pi__r

# Ideal gas part
phi0 = ln(delta)+_n0[1]+_n0[2]*tau+_n0[3]*ln(tau)
for i in range(4,9):
    phi0 += _n0[i]*ln(1-exp(-_gamma0[i]*tau))

Psi={}; Theta={}; Delta={}
for i in range(55,57):
    Psi[i]   = exp(-_C[i]*(delta-1)**2-_D[i]*(tau-1)**2)
    Theta[i] = (1-tau)+_A[i]*((delta-1)**2)**(1/(2*_beta[i]))
    Delta[i] = Theta[i]**2+_B[i]*((delta-1)**2)**_a[i]


# Eq 6.6 p. 429
# Residual part
phi__r = 0
for i in range(1,8):
    phi__r += _n[i]*delta**_d[i]*tau**_t[i]
for i in range(8,52):
    phi__r += _n[i]*delta**_d[i]*tau**_t[i]*exp(-delta**_c[i])
for i in range(52,55):
    phi__r += _n[i]*delta**_d[i]*tau**_t[i]*exp(-_alpha[i]*(delta-_epsilon[i])**2-_beta[i]*(tau-_gamma[i])**2)
for i in range(55,57):
    phi__r += _n[i]*Delta[i]**_b[i]*delta*Psi[i]

dphi__rdrho = diff(phi__r,rho)
dphi__rddelta = rho_c * dphi__rdrho
pressure_relation = P/(rho*R*T)-1-delta*dphi__rddelta

# <codecell>


def get_water_density(P_val=None, T_val=None):
    if not P_val: P_val = sympify(101.3e3) * units.pascal
    if not T_val: T_val = sympify(298.15)  * units.kelvin

    # If P is without unit:
    try:
        float(P_val/units.pascal)
    except:
        P_val *= units.pascal

    # If T is without unit
    try:
        float(T_val/units.kelvin)
    except:
        T_val *= units.kelvin

    def f0(x_rho):
        return pressure_relation.subs({P:P_val,T:T_val,rho:x_rho}).evalf()

    rho0=sympify(1000.00)
    return find_root(f0,rho0,unit=units.kg/units.meter**3,maxiter=25,verbose=False)

#print a.evalf()
#solve(pressure_relation, rho)

P0=sympify(101.3e3)*units.pascal
Tdash=sympify(298.15)*units.kelvin

get_water_density(P0,Tdash)

rho0 = sympify(1000.0)*units.kg/units.meter**3
pressure_relation.subs({P:P0,T:Tdash,rho:rho0}).evalf()


# For verification:
# From table 6 in IAPWS-Rev.pdf
abstol    = sympify(1e-7)
T_ver     = sympify(500.0) * units.kelvin
rho_ver   = sympify(838.025) * units.kg / units.meter**3
subs_dict = {T:T_ver,rho:rho_ver}



dphi0ddelta       = diff(phi0,rho)     *rho_c
dphi__rddelta     = diff(phi__r,rho)   *rho_c
dphi0dtau         = diff(phi0,T)       *T_c
dphi__rdtau       = diff(phi__r,T)     *T_c
d2phi0ddelta2     = diff(phi0,rho,2)   *rho_c**2
d2phi__rddelta2   = diff(phi__r,rho,2) *rho_c**2
d2phi0dtau2       = diff(phi0,T,2)     *T_c**2
d2phi__rdtau2     = diff(phi__r,T,2)   *T_c**2
d2phi0ddeltatau   = diff(phi0,rho,T)   *T_c*rho_c
d2phi__rddeltatau = diff(phi__r,rho,T) *T_c*rho_c


# assert(abs(            phi0.subs() - phi0_val              ) < abstol)
# assert(abs(          phi__r.subs({T:T_val,rho:rho_val}) - phi__r_val            ) < abstol)
# assert(abs(     dphi0ddelta.subs({T:T_val,rho:rho_val}) - phi0_delta_val        ) < abstol)
# assert(abs(   dphi__rddelta.subs({T:T_val,rho:rho_val}) - phi__r_delta_val      ) < abstol)
# assert(abs(       dphi0dtau.subs({T:T_val,rho:rho_val}) - phi0_tau_val          ) < abstol)
# assert(abs(     dphi__rdtau.subs({T:T_val,rho:rho_val}) - phi__r_tau_val        ) < abstol)
# assert(abs(    dphi0ddelta2.subs({T:T_val,rho:rho_val}) - phi0_deltadelta_val   ) < abstol)
# assert(abs(  dphi__rddelta2.subs({T:T_val,rho:rho_val}) - phi__r_deltadelta_val ) < abstol)
# assert(abs(      dphi0dtau2.subs({T:T_val,rho:rho_val}) - phi0_tautau_val       ) < abstol)
# assert(abs(    dphi__rdtau2.subs({T:T_val,rho:rho_val}) - phi__r_tautau_val     ) < abstol)
# assert(abs(  dphi0ddeltatau.subs({T:T_val,rho:rho_val}) - phi0_deltatau_val     ) < abstol)
# assert(abs(dphi__rddeltatau.subs({T:T_val,rho:rho_val}) - phi__r_deltatau_val   ) <       abstol)


# print dphi0dtau.subs({T:T_val,rho:rho_val})
# print phi0_tau_val


# Test triple point
T_t=sympify(273.16)*units.kelvin



def staurated_liquid_density():
    raise NotImplemented # Not finished..
    # Saturated liquid density
    rho__prime = T*dp_sigmadT/Beta # Eq 2.3 p 398
    v = (1-T/T_c)
    b1 =  1.99274064
    b2 =  1.09965342
    b3 = -0.510839303
    b4 = -1.75493479
    b5 = -45.5170352
    b6 = -6.74694450e5
    delta_prime = 1 + b1*v**(1.0/3.0) + \
                          b2*v**(2.0/3.0) + \
                          b3*v**(5.0/3.0) + \
                          b4*v**(16.0/3.0) + \
                          b5*v**(43.0/3.0) + \
                          b6*v**(110.0/3.0)

