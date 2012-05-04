
def get_unit(arg):
    from sympy.physics import units
    if hasattr(arg, 'as_ordered_factors'):
	found_units = []
	for factor in arg.as_ordered_factors():
	    if str(factor) in units.__dict__:
		found_units.append(factor)
	    else:
		if hasattr(factor, 'base'):
		    if str(factor.base) in units.__dict__:
			found_units.append(factor)
	return reduce(mul,found_units)
    else:
	return 1.0

def get_unitless(arg):
    return arg / get_unit(arg)
