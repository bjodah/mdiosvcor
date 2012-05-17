from sympy import sympify
from operator import mul
from functools import reduce, wraps

# TODO: improve pickle_cached memoization to support critical and non_critial arguments, and support *args, **kwargs combinations.

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
	if len(found_units) > 0:
	    return reduce(mul,found_units)
	else:
	    return 1.0
    else:
	return 1.0

def get_unitless(arg):
    return arg / get_unit(arg)


def get_sympified(instance):
    if hasattr(instance, 'iteritems'):
	return dict([(k,sympify(v)) for k,v in instance.iteritems()])
    elif hasattr(instance, '__iter__'):
	return [sympify(x) for x in instance]
    else:
	NotImplemented


def memoize(f):
    """
    Decorator to enable caching for computationally heavy functions
    """
    cache={}
    @wraps(f)
    def helper(*args):
        if args not in cache:
            cache[args]=f(*args)
        return cache[args]
    return helper


def pickle_cached(f):
    import cPickle as pickle
    import os.path
    from os import environ
    dirpath = environ.get('PICKLE_CACHE_DIR','.')
    fname = '.pickle_cached__'+f.__name__
    path = os.path.join(dirpath, fname)
    if os.path.exists(path):
        cache = pickle.load(open(path,'rb'))
    else:
        cache = {}
    @wraps(f)
    def helper(*args):
        if args not in cache:
            cache[args]=f(*args)
            pickle.dump(cache, open(path, 'wb'))
        return cache[args]
    return helper


class ParameterStore(object):
    """
    A convenience class for storing parameters with
    values and then at later point in time choose
    whether or not the retrieve them unitless.
    """

    def __init__(self, parameters, return_unitless=False):
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

