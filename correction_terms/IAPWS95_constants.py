#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Provides constants used in the IAPWS95 formulation.

# Author: Björn Dahlgren
# Implemented for use in research project in the IGC group at ETH

# To the extent possible under law, Bjoern Dahlgren has waived all
# copyright and related or neighboring rights to this work.

from sympy_helpers import get_unit, get_unitless
from sympy import sympify
from sympy.physics import units
# Constants used in the master equation
class ParameterStore(object):
    """

    """

    def __init__(self, parameters, return_unitless=True):
        """

        Arguments:
        - `unitless`: Default
        """
        self._parameters = parameters
        self.return_unitless = return_unitless

    def keys(self): return self._parameters.keys()

    def __getitem__(self, item):
        if self.return_unitless:
            return get_unitless(self._parameters[item])
        else:
            return self._parameters[item]

ref_variables = ParameterStore({'T_c': sympify(647.096) * units.kelvin,
                                     'rho_c': sympify(322.0) * units.kg / units.meter**3,
                                     'R': sympify(461.51805) * units.joule / units.kelvin / units.kg})

standard_state_variables = ParameterStore({'P0'    : sympify(101.325e3)*units.pascal,
                                           'Tdash' :	sympify(298.15)*units.kelvin
                                           })


_n0 = [0.0,
       -8.3204464837497, 6.6832105275932, 3.00632,
        0.012436,        0.97315,         1.27950,
        0.96956,         0.24873]


_gamma0 = [0.0,
           0.0,               0.0,               0.0,
           1.28728967,      3.53734222,      7.74073708,
           9.24437796,     27.5075105]


# The raw data is copied from IAPWS95.py by Kiran Pashikanti
_n = [0.0,
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
         -0.14874640856724e+0,  0.31806110878444e+0]

_c = [0,
          0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
          1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
          2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
          2, 2, 2, 3, 3, 3, 3, 4, 6, 6, 6, 6]
_d =  [0,
          1, 1, 1, 2,  2,  3,  4,  1,  1, 1, 2,  2,  3,  4,
          4, 5, 7, 9, 10, 11, 13, 15,  1, 2, 2,  2,  3,  4,
          4, 4, 5, 6,  6,  7,  9,  9,  9, 9, 9, 10, 10, 12,
          3, 4, 4, 5, 14,  3,  6,  6,  6, 3, 3,  3]
_t =  [0,
          -0.5, 0.875,  1.0,  0.5,  0.75, 0.375,  1.0,  4.0,  6.0, 12.0,  1.0,
           5.0, 4.0  ,  2.0, 13.0,  9.0 , 3.0  ,  4.0, 11.0,  4.0, 13.0,  1.0,
           7.0, 1.0  ,  9.0, 10.0, 10.0 , 3.0  ,  7.0, 10.0, 10.0,  6.0, 10.0,
          10.0, 1.0  ,  2.0,  3.0,  4.0 , 8.0  ,  6.0,  9.0,  8.0, 16.0, 22.0,
          23.0,23.0  , 10.0, 50.0, 44.0, 46.0  , 50.0,  0.0,  1.0,  4.0]

#offset constants, no point in creating giant nearly empty lists
_alpha	 = {52: 20.0,  53:  20.0,  54:  20.0}
_beta	 = {52: 150.0, 53: 150.0,  54: 250.0, 55:0.3, 56:0.3}
_gamma	 = {52: 1.21,  53:   1.21, 54:   1.25}
_epsilon = {52: 1.0,   53:   1.0,  54:   1.0 }

_a		 = {55:   3.5,   56:   3.5}
_b		 = {55:   0.85,  56:   0.95}
_B		 = {55:   0.2,   56:   0.2}
_C		 = {55:  28.0,   56:  32.0}
_D		 = {55: 700.0,   56: 800.0}
_A		 = {55:   0.32,  56:   0.32}

