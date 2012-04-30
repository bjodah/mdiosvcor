#!/usr/bin/env python
# -*- coding: utf-8 -*-

from manuscript_equations import get_Delta_G_LS

P_val = 5000   * units.bar
T_val = 310.65 * units.kelvin
N_W   = 1024          # Number of water molecules in MD simulation
ion   = 'sod'         # What ion is used (sod: Na^+ / cls: Cl^-)

print get_Delta_G_LS(P_val, T_val, N_W, ion, 'all')
