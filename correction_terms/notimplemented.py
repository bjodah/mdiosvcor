# Entities below are currently not used. But are kept for possible
# future reuse of the code and constants
# ==================================================================
# ==================================================================



# Other constants never used explicitly
# =====================================
# Critical pressure
# P_c = 22.064e6 * units.pascal

# Triple point,
T_t	    = sympify(273.16)*units.kelvin    # Eq 2.1a p. 398
p_t	    = sympify(611.657)*units.pascal   # Eq 2.1b p. 398
rho_prime_t = sympify(999.793)*units.kg/units.meter**3
rho_bis_t   = sympify(0.00485458) * units.kg / units.meter**3



def staurated_liquid_density():
    raise NotImplemented # Not finished..
    # Saturated liquid density
    rho__prime = T*dp_sigmadT/Beta # Eq 2.3 p 398
    v = (1-T/T_c)
    b  =  [0,
           1.99274064,
           1.09965342,
          -0.510839303,
          -1.75493479,
         -45.5170352,
          -6.74694450e5]
    delta_prime = 1 + b[1]*v**(1.0/3.0) + \
                      b[2]*v**(2.0/3.0) + \
                      b[3]*v**(5.0/3.0) + \
                      b[4]*v**(16.0/3.0) + \
                      b[5]*v**(43.0/3.0) + \
                      b[6]*v**(110.0/3.0)

