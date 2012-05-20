from sympy import sympify
from operator import mul
from functools import reduce, wraps

# TODO: improve pickle_cached memoization to support critical and non_critial arguments, and support *args, **kwargs combinations.

def get_unit(arg):
    """
    Returns what units (in sympy.physics.units)
    an expression has (terms of mixed units not supported).
    A unitless expression returns 1.0 as unit.
    """
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
    """
    A more advanced memoization decorator which stores
    calculated values to file. Both the return value of
    the function and the output to STDOUT is cached.
    """
    import cPickle as pickle
    import os.path
    from os import environ
    from collections import defaultdict
    dirpath = environ.get('PICKLE_CACHE_DIR', '.')
    fname = '.pickle_cached__'+f.__name__
    path = os.path.join(dirpath, fname)
    if os.path.exists(path):
        cache = pickle.load(open(path, 'rb'))
    else:
        cache = defaultdict(dict)
    @wraps(f)
    def wrapper(*args):
        if args not in cache:
	    OutputCache = TeeOutputCacher()
            cache[args]['retval']=f(*args)
	    cache[args]['stdout']=OutputCache.content
	    del OutputCache
            pickle.dump(cache, open(path, 'wb'))
	else:
	    print ''.join(cache[args]['stdout']),
        return cache[args]['retval']
    return wrapper


import sys

class OutputCacher(object):
    """
    By initializing a OutputCacher instance
    sys.stdout is cached in that instance until it
    is deleted. (Remember to dump content before deletion).
    Only one instance is allowed to catch unexpected behaviour
    """
    # http://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-stdout-to-a-log-file-in-python
    # http://stefaanlippens.net/redirect_python_print
    num_instances = 0
    def __init__(self):
        self.content = []
        self.stdout  = sys.stdout
        sys.stdout   = self
	if OutputCacher.num_instances == 0:
	    OutputCacher.num_instances = 1
	else:
	    raise ValueError('Tried to initalize multiple OutputCachers')

    def __del__(self):
        sys.stdout = self.stdout
	OutputCacher.num_instances -= 1

    def write(self, data):
        self.content.append(data)

    def dump(self):
        return ''.join(self.content)


class TeeOutputCacher(OutputCacher):
    """
    Extends OutputCacher and instead of rediriecting
    sys.stdout it is copied to cache.
    (useful for long running functions which prints
    information over time)
    """
    def write(self, data):
        self.content.append(data)
	self.stdout.write(self.dump())


class ParameterStore(object):
    """
    A convenience class for storing parameters with
    values (with units) and then at later point in time choose
    to retrieve them with or without units.
    """

    def __init__(self, parameters, return_unitless=False):
        """
        Arguments:
        - `parameters`:      A dict of paramaters
	- `return_unitless`: Default: False
        """
        self._parameters = parameters
        self.return_unitless = return_unitless

    def keys(self): return self._parameters.keys()

    def __getitem__(self, item):
        if self.return_unitless:
            return get_unitless(self._parameters[item])
        else:
            return self._parameters[item]
