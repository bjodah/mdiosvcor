#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This work is open source and is released under the
# 2-clause BSD license (see LICENSE.txt for further information)
# Copyright (c) 2011, 2012, Björn Ingvar Dahlgren

import sys
import cPickle as pickle
import os

from sympy import sympify
from operator import mul
from functools import reduce, wraps
from tempfile import gettempdir

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

def in_terms_of(self, terms):
    """
    From http://stackoverflow.com/questions/2038100/sympy-how-to-return-an-expression-in-terms-of-other-expressions

    It would be soo much nicer if sympy would support
    the PyPI package `Quantities`

    Exapmle usage:

        >>> x = in_terms_of('J',['gram','mile','hour'])
        >>> x
        '9765625000*mile**2*gram/(1951609*hour**2)'
        >>> in_terms_of(x,['J'])
        'J'

    """
    from sympy.physics import units

    expr2 = eval(str(self), units.__dict__)
    converter = {}
    for term in terms:
        term_expr = units.__dict__[term]
	print term_expr
        coeff = term_expr.as_coeff_terms()[0]
        unit = units.Unit('my_'+term, term)
        converter[term_expr/coeff] = unit/coeff
    return str(expr2.subs(converter))



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
    def wrapper(*args):
        if args not in cache:
            cache[args]=f(*args)
        return cache[args]
    return wrapper


class PickleDict(dict):
    """
    A dictionary with attribute-style access. It maps attribute access to
    the real dictionary.

    The cache file must be set before invoking `dump`
    (or assigning key/value pairs to the instance if autodump is enabled).
    It can be done after initialization by invoking `set_cache_file`

    The code is largely based on:
    http://stackoverflow.com/questions/9409045/
         how-to-create-a-persistant-class-using-pickle-in-python

    TODO: Improve only to write changes to file instead of redumping entire dictionary.
    """

    def __init__(self, *args, **kwargs):
        self._autodump = False # In order to dump cache file needs to be initialized
        super(self.__class__, self).__init__(self, *args, **kwargs)

    def set_cache_file(self, cache_name, cache_dir_path=None,
                            fname_prefix=None, autoload=True,
                            autodump=True):
        """
        Sets the file path of the pickle file.
        If `cache_dir_path` is not specified, tempfile.gettempdir will set one.

        By leaving the kwarg `autoload` at its default (True) will load the
        pickle file if it exists
        """
        self._autodump = autodump
        if fname_prefix == None:
            fname_prefix = '.'+self.__class__.__name__+'_'
        if not cache_dir_path:
            cache_dir_path = gettempdir()
        fname = fname_prefix+cache_name
        self._cache_file = os.path.join(cache_dir_path, fname)
        if autoload:
            if os.path.exists(self._cache_file):
                self.load()


    def dump(self):
        pickle.dump(self, open(self._cache_file, 'wb'))

    def load(self):
        self.__setstate__(pickle.load(open(self._cache_file, 'rb')).__getstate__())

    def unlink(self):
        if os.path.exists(self._cache_file): os.unlink(self._cache_file)

    def __getstate__(self):
        # Only pickle the underlying dict, not the PickleCache instance
        return self.items()

    def __setstate__(self, items):
        # When unpickling, we can assume reading an ordinary dictionary
        # The meta-data which makes the PicklCache instance is already present
        for key, val in items:
            self[key] = val

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, dict.__repr__(self))

    def __setitem__(self, key, value):
        # Not idempotent, dump data if autodump enabled
        retval = super(self.__class__, self).__setitem__(key, value)
        if hasattr(self, '_autodump'):
            if self._autodump: self.dump()
        return retval

    def __getitem__(self, name):
        # Idempotent operation, don't dump even if autodump is enabled
        return super(self.__class__, self).__getitem__(name)

    def __delitem__(self, name):
        # Not idempotent, dump data if autodump enabled
        retval = super(self.__class__, self).__delitem__(name)
        if self._autodump: self.dump()
        return retval

    def copy(self):
        return self.__class__(self)

def test_PickleDict():
    # Tests the functionality of the PickleDict class
    pd = PickleDict()
    pd.set_cache_file(cache_name='test', cache_dir_path='/tmp/',
                           autoload=False)
    pd['foo'] = ('a', 42)
    pd['bar'] = ('b',137)
    print pd['foo']
    print pd.get('bar','oops no bar!')
    print pd.get('baz','oops no baz!')
    del pd

    pd = PickleDict()
    pd.set_cache_file(cache_name='test',cache_dir_path='/tmp/')
    assert dict(pd) == {'foo': ('a',42), 'bar': ('b',137)}
    pd.unlink()
    del pd



def adv_memoize(cache_name=None, cache_stdout=True, cache_dir_path=None, autoload=True, verbose=False):
    """
    A more advanced memoization decorator factory which
    creates decorators that enable memoization which stores
    calculated values to file. Both the return value of
    the function and, optionally, the output to STDOUT is cached.

    If not the cache_dir_path kwarg is given the environment variable
    `MEMOIZE_CAHCE_DIR` will be looked for (fallback is then PWD).

    TODO: Use python standard library logging module instead
    TODO: Make it more compatible with functools.lru_cache from Python 3.2+
    """
    if cache_dir_path == None:
        cache_dir_path = os.environ.get('MEMOIZE_CACHE_DIR', '.')
    def decorator(f):
        if cache_name == None:
            _cache_name = f.__name__
        else:
            _cache_name = cache_name

        # Set keyword arguments to `set_cache_file` CacheDict instance method
        scf_kwargs={'cache_dir_path' : cache_dir_path,
                    'fname_prefix'   : ".adv_memoize_cache_",
                    'autoload'       : autoload,
                    'autodump'       : True}

        retval_cache = PickleDict()
        retval_cache.set_cache_file(f.__name__+'_retval',**scf_kwargs)

        if cache_stdout:
            output_cache = PickleDict()
            output_cache.set_cache_file(f.__name__+'_output',**scf_kwargs)

        @wraps(f)
        def wrapper(*args):
            if args not in retval_cache:
                # Function output not in cache
                if cache_stdout:
                    # Let's copy stdout to cache
                    with TeeOutputCacher() as cacher:
                        retval_cache[args] = f(*args)
                        output_cache[args] = cacher.dump()
                else:
                    retval_cache[args] = f(*args)
            else:
		if verbose: print "Cached value found"
                if cache_stdout:
                    # Print cached output to stdout
                    print output_cache[args],
            return retval_cache[args]
        return wrapper
    return decorator

def test_adv_memoize__1():

    @adv_memoize(cache_stdout=False, autoload=False)
    def square(x):

        print "Calculating square of " + str(x)
       # print  >> sys.stderr, "DEBUGINFO: x=" + str(x)
        square.i += 1
        return x**2

    square.i = 0
    square(2)
    square(2)
    square(2)
    assert square.i == 1

def test_adv_memoize__2():
    @adv_memoize(cache_stdout=True, autoload=False)
    def cube(x):
        print "Calculating cube of " + str(x)
        cube.i += 1
        return x**3

    cube.i = 0
    cube(2)
    cube(2)
    cube(2)
    assert cube.i == 1

def test_adv_memoize__3():

    @adv_memoize(cache_stdout=True, autoload=False)
    def root(x):

        print "Calculating root of " + str(x)
        print  >> sys.stderr, "DEBUGINFO: x=" + str(x)
        root.i += 1
        return x**0.5

    root.i = 0
    root(2)
    root(2)
    root(2)
    assert root.i == 1



class OutputCacher(object):
    """
    By using a OutputCacher instance in a `with` statement
    sys.stdout is cached in that instance until the with statement
    is completed. (Remember to dump content before exiting with statement).
    """
    # http://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-stdout-to-a-log-file-in-python
    # http://stefaanlippens.net/redirect_python_print

    def __init__(self):
        self.content = []

    def __enter__(self):
        # Used in `with` statement
        self.stdout  = sys.stdout
        sys.stdout   = self
        return self

    def __exit__(self, type, value, traceback):
        sys.stdout = self.stdout

    def write(self, data):
        self.content.append(data)

    def dump(self):
        return ''.join(self.content)

def test_OutputCacher():
    print 'Creating OutputCacher instance'
    with OutputCacher() as oc:
        print 'hello'
        print 'world'
        tmp = oc.dump()

    print 'Just exited with statement, lets look what was put into output:'
    print tmp,
    print 'Well, I hope it reads \"hello\nworld\"'


class TeeOutputCacher(OutputCacher):
    """
    Extends OutputCacher and instead of rediriecting
    sys.stdout it is copied to cache.
    (useful for long running functions which prints
    information over time)
    """
    def write(self, data):
        self.content.append(data)
	self.stdout.write(data)


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

def get_dict_combinations_for_diff(d):
    """
    Useful for generating a list ordered for differentiating
    a function step by step to higher and higher orders with
    respect to a variable.

    E.g.

    In [14]: combo.get_dict_combinations_for_diff({'x':3,'y':2})
    Out[14]:
    [{'x': 0, 'y': 0},
     {'x': 0, 'y': 1},
     {'x': 0, 'y': 2},
     {'x': 1, 'y': 0},
     {'x': 1, 'y': 1},
     {'x': 1, 'y': 2},
     {'x': 2, 'y': 0},
     {'x': 2, 'y': 1},
     {'x': 2, 'y': 2},
     {'x': 3, 'y': 0},
     {'x': 3, 'y': 1},
     {'x': 3, 'y': 2}]

    """

    from itertools import product

    list_of_dicts = []
    ranges = []
    keys = d.keys()
    for key in keys:
        ranges.append(range(0,d[key]+1))

    for cmb in product(*ranges):
        list_of_dicts.append(dict(zip(keys, cmb)))

    return list_of_dicts


