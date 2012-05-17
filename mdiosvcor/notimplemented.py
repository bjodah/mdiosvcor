
# Entities below are currently not used. But are kept for possible
# future reuse of the code and constants. This is more or less a
# trash can.
# ==================================================================
# ==================================================================


# Code for handling unitless evaluation of taylor expansion

TayParams = namedtuple('TayParams', ['f0', 'dfdP', 'dfdT', 'd2fdP2', 'd2fdPdT', 'd2fdT2'])

def TayExpr(x, y, f0_val, dfdx_val, dfdy_val,
	    d2fdx2_val, d2fdxdy_val, d2fdy2_val,
	    x0=0.0, y0=0.0):
    """
    Returns a 2nd order Taylor expansion expression
    of a single valued function of two variables (x,y).
    It is built from derivative values and x0.
    """
    f0,dfdx,dfdy,d2fdx2,d2fdxdy,d2fdy2 = \
		  symbols('f0 dfdx dfdy d2fdx2 d2fdxdy d2fdy2')
    return f0 + dfdx*(x-x0) + dfdy*(y-y0) + \
       1/2*(d2fdx2*(x-x0)**2 + \
	       2*d2fdxdy*(x-x0)*(y-y0) + \
	       d2fdy2*(y-y0)**2)

class TayExprFactory(object):
    """

    """

    def __init__(self, f0_val, dfdx_val, dfdy_val,
	    d2fdx2_val, d2fdxdy_val, d2fdy2_val,
	    x0=0.0, y0=0.0, return_unitless=False):
        """

        """
        self._parameters = {
	    'f0_val'	  : f0_val,
	    'dfdx_val'	  : dfdx_val,
	    'dfdy_val'	  : dfdy_val,
	    'd2fdx2_val'  : d2fdx2_val,
	    'd2fdxdy_val' : d2fdxdy_val,
	    'd2fdy2_val'  : d2fdy2_val,
	    'x0'	  : x0,
	    'y0'	  : y0,
	    }
        self.return_unitless = return_unitless

    def keys(self): return self._parameters.keys()

    def __getitem__(self, item):
        if self.return_unitless:
            return get_unitless(self._parameters[item])
        else:
            return self._parameters[item]

    def __call__(self, return_unitless=None):
	if not return_unitleess:
	    return_unitless = self.return_unitless
	return TayExpr(**self)





# IAPWS 95 specific
############################################################

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

