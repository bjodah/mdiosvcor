#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This work is open source and is released under the
# 2-clause BSD license (see LICENSE.txt for further information)
# Copyright (c) 2011, 2012, Björn Ingvar Dahlgren

"""
Module implementing common root finding algorithms.
Implemented with sympy and sympy.physics.units compability in mind.
"""

from functools import reduce
from operator import and_, mul
from prj_helpers import get_unit, get_unitless
from sympy.utilities.autowrap import autowrap

try:
    import numexpr as ne
except ImportError:
    ne = False


def find_root(func,
              x0,
              preferred_method=None,
              dfdx=None,
              dx0=None,
              xunit=1.0,
	      yunit=1.0,
              abstol=None,
              yabstol=None,
              xabstol=None,
	      xreltol=None,
              maxiter=25,
              return_intermediate_results=False,
              verbose=False):
    """
    If no preferred method is specified:
        If dfdx is specified, assume that Newtons method is wanted.
        Else assume secant method.

    TODO: let unit be automatically taken from x0 and func
    """

    if preferred_method == "newton" or \
        (not preferred_method and dfdx):
        if verbose: print "Using Newton's method."
        method_cb, method_args = (newton_generator, (func, x0, dfdx))
    elif preferred_method == "secant" or \
        (not preferred_method and not dfdx):
        if verbose: print "Using the secant method."
        if not dx0:
            # If no dx0 specified take an arbitrary step forward
	    if x0 != 0:
		dx0 = x0*1e-1
	    else:
		dx0 = 1e-7*xunit
            if verbose: print "Assuming dx0 = {0}".format(dx0)
        method_cb, method_args = (secant_generator, (func, x0, dx0))
    else:
        NotImplemented

    # Setup tolerance requierments
    satisfication_functions = []

    def satisfied():
        satisfication = [sf() for sf in satisfication_functions]
        return reduce(and_,satisfication)

    def satisfied_yabstol():
        return abs(dy/yunit)<yabstol

    def satisfied_xabstol():
        return abs(dx/xunit)<xabstol

    def satisfied_xreltol():
        return abs(dx/x)<xreltol


    if abstol:
        if xabstol or yabstol:
            raise ValueError("Both abstol and xabstol and/or yabstol specified")
        xabstol = abstol
        yabstol = abstol


    if not yabstol and not xabstol and not xreltol:
        xreltol = 1e-8
        if verbose: print "Since no tolerances were specified," + \
           " xreltol assumed = {0:12.5g}".format(xreltol)

    if yabstol: satisfication_functions.append(satisfied_yabstol)
    if xabstol: satisfication_functions.append(satisfied_xabstol)
    if xreltol: satisfication_functions.append(satisfied_xreltol)

    # Handle itermediate results if requested
    if return_intermediate_results:
        interm_res = []

    def append_res(x,y):
	interm_res.append((x,y))


    # Initialize the generator of chosen method
    method_gen = method_cb(*method_args)
    i = 0
    x0, y0 = method_gen.next()
    dx, dx_old, dy = None, None, None
    if verbose: print_info(i,x0,dx,y0,dy,dx_old,xunit,yunit)
    for x,y in method_gen:
        i     += 1
        dx_old = dx
	dy_old = dy
        dx, dy = x - x0, y - y0
        if verbose: print_info(i,x,dx,y,dy,dx_old,xunit,yunit)
        if return_intermediate_results: append_res(x,y)
        if satisfied(): break
        if i == maxiter:
            raise ValueError("Maximum number of iterations reached")
        x0, y0 = x, y

    if return_intermediate_results:
        return zip(*interm_res)
    else:
	return x, dx


def secant_generator(f, x0, dx0):
    """
    Recursion formula for the "Secant method"
    """
    x1, x2 = x0+dx0, x0
    f1, f2 = f(x1), f(x2)
    yield x2, f2
    yield x1, f1
    while True:
        x = x1 - f1*(x1-x2)/(f1-f2)
        x1, x2 = x, x1
        f2, f1 = f1, f(x1)
        yield x1, f1


def newton_generator(f, x0, dfdx):
    """
    Newton's method or more correctly: Newton-Raphson's method
    """
    while True:
        dfdx0 = dfdx(x0)
        f0 = f(x0)
        yield x0, f0
        x0 -= f0/dfdx0


def fixed_point_iteration_generator(G, x0, dGdx=None):
    """
    Usually not the best option but if a fixed point iteration
    is needed this generator provides the recursion formula.
    """
    conv_limit = 1 - 1e-6
    if dGdx:
        dGdx0 = dGdx(x0)
        if abs(dGdx0) > conv_limit:
            raise ValueError("This equation will not converge")
        elif abs(dGdx0) > 0.5:
            print "Warning: equation probably not optimal for fix point iteration"
    G_old = x0
    G_new = G(G_old)
    dx_old = G_new - G_old
    yield G_new
    while True:
        G_new = G(x)
        dx = G_new - G_old
        m = dx / dx_old
        if abs(m) > conv_limit:
            raise ValueError("This equation will not converge")
        yield G_new
        G_old = G0
        dx_old = dx


def print_info(i,x,dx,y,dy,dx_old, xunit=1.0, yunit=1.0,
                 fmtstr="{0: >3}{1: >27}{2: >27}"+\
                        "{3: >27}{4: >27}{5: >27}"):
    """
    For studying the convergence it is helpful to
    see how the values changes between steps.
    Call find_root with "verbose=True" in order
    to enable printing of intermediate values.
    """
    if xunit != 1.0 or yunit != 1.0:
	from sympy.mpmath import mp
	mp.dps = 10

    if i==0:
	if xunit != 1.0:
	    sx	      = "x [{0}]".format(xunit)
	    sdx	      = "dx [{0}]".format(xunit)
	    sdxidxim1 = "dxi/dx(i-1)**2 [{0}]".format(1/xunit)
	else:
	    sx	      = "x"
	    sdx	      = "dx"
	    sdxidxim1 = "dxi/dx(i-1)**2"
	if yunit != 1.0:
	    sy	      = "y [{0}]".format(yunit)
	    sdy	      = "dy [{0}]".format(yunit)
	else:
	    sy	      = 'y'
	    sdy	      = 'dy'

        print fmtstr.format("i",sx,sdx,sy,sdy,sdxidxim1)
        print fmtstr.format(i, x/xunit, "-", y/yunit, "-", "-")
    elif i ==1:
        print fmtstr.format(i, x/xunit, dx/xunit, y/yunit, dy/yunit, '-')
    else:
        print fmtstr.format(i, x/xunit, dx/xunit, y/yunit, dy/yunit,
			    dx/dx_old**2*xunit)


def solve_relation_num(rel,
		       subsd,
		       varied,
		       initial_guess,
		       use_numexpr=False, # Info to get_cb_from_rel
		       unitless=False,    # Info to get_cb_from_rel
		       bincompile=False,  # Info to get_cb_from_rel
		       verbose=False,
		       **kwargs # Additional options to find_root
		       ):
    """
    Solves non-linear (sympy) equation numerically
    """
    f0 = get_cb_from_rel(rel, subsd, varied, use_numexpr,
			 unitless, bincompile)
    return find_root(f0, initial_guess, verbose=verbose, **kwargs)


def test_solve_relation_num(use_numexpr=False, bincompile=False):
    from sympy import symbols, Function, ln
    x,y = symbols('x,y')
    f = symbols('f',cls=Function)(x,y)
    relation = ln(x*y+f)-f
    subsd = {x:2,y:1}
    varied = f
    initial_guess = 1.0
    unitless = True
    x, y = solve_relation_num(relation, subsd, varied,
			     initial_guess, use_numexpr,
			     unitless, bincompile,
			     verbose=True, yabstol=1e-9)

    assert abs(x-1.14619322062) < 1e-8

    x,y = symbols('x,y')
    g = symbols('g', cls=Function)(x,y)
    relation = x**2+y**2-g # Circle
    g_val,delta_rel = solve_relation_num(relation, {x:1, y:0}, g, 0.1, use_numexpr, verbose=True, abstol=1e-8)
    assert abs(g_val-1) < 1e-7


def solve_relation_for_derivatives(rel,
				   subsd,
				   func,  # Function which to solve for (incl. derivs)
				   initial_guess_vals,
				   diff_wrt={},
				   diff_wrt_map={},
				   use_numexpr=False,  # Info to get_cb_from_rel
				   unitless=False,     # Info to get_cb_from_rel
				   bincompile=False,   # Info to get_cb_from_rel
				   verbose=False,
				   **kwargs
				   ):
    """

    E.g. assume we are trying to find the value of:
     f, dfdx, d2fdx2, dfdy, d2fdy2 for the relation:

     f(x,y) = ln(x*y+f(x,y))

     f = symbols('f', cls=Function)(x,y)
     rel = ln(x*y+f) - f   (x*y >= 1)

     variables = (x,y)

    Ideas for improvements: from iterations use finite difference to
    generate a better starting guess for next derivative

    """

    from combo import get_dict_combinations_for_diff
    from operator import add
    from functools import reduce
    from sympy import Derivative

    drel      = {} # differentiated relaion
    deriv     = {} # symbolic derivative sought for
    deriv_val = {} # Value of the derivative sought for
    err_deriv_val = {}

    # Start by differentiating with respect to one variable to first order,
    # then succesively increase order, then switch variable.
    for diff_step in get_dict_combinations_for_diff(diff_wrt):
	if verbose: print "Solving for derivative with respect to: "+\
	   str(diff_step)
        # Lets make a signature e.g. ((x,1), (y,1)) corresponds to d2fdxdy
        signature = tuple(diff_step.items()) # We need something hashable

        # Now let us differentiate the relation wrt the variables of the sig.
        drel[signature] = rel.diff(*reduce(add,signature))

        # Let us identify the symbolic form of the sought for derivative
        deriv[signature] = Derivative(func, *reduce(add,signature))


        # Find the relation by feeding the derived relation to "solve_relation_num"
        # using the symbolic derivative deriv[signature] as unknown.
        # But first we need a starting guess.
        #
	if signature in initial_guess_vals:
            initial_guess = initial_guess_vals[signature]
	    if verbose: print "Using provided inital guess: ", str(initial_guess)
	    kwargs.update({'xunit':get_unit(initial_guess)})
        else:
            # There exist a 'parent' derivative from which we can
            # extract an initial guess for "solve_relation_num"
	    if verbose: print "No initial guess provided for this derivative, using forward difference estimate from one step order derivative: "
            parent_sig = None
            for k,v in diff_step.iteritems():
                if v > 0:
                    mod_diff_step = diff_step
                    mod_diff_step[k] = v - 1
                    if tuple(mod_diff_step.items()) in deriv_val.keys():
                        parent_sig = tuple(mod_diff_step.items())
                        break
            if parent_sig == None: raise ValueError('no parent signature found!')

            # Now let us identify which variable which derivation order
            # has increased by 1 from parentrelation.

            new_order_var = None
            for candidate in diff_step.keys():
		candidate_dict = dict(parent_sig)
		candidate_dict[candidate] += 1
                if candidate_dict == dict(signature):
		    new_order_var = candidate#diff_wrt_map[candidate]
	    if new_order_var == None: raise ValueError('The varible was not correctly determined')
	    if get_unitless(subsd[new_order_var]) > 1e-12:
		h = subsd[new_order_var]*1e-1
	    else:
		h = 1e-6 # Arbitrary step
		if not unitless: h *= get_unit(subsd[new_order_var])
	    if 'dx0' in kwargs: kwargs.pop('dx0')
	    subsd_ph = dict(subsd.items()) # Deepcopy imitation
	    subsd_ph[new_order_var] += h
	    # Forward difference
            initial_guess = (drel[parent_sig].subs(subsd_ph) - \
			     drel[parent_sig].subs(subsd)) / h#*get_unit(initial_guess_func_vals[])
	    kwargs.update({'xunit':get_unit(initial_guess)})

	if verbose: print 'Determining:', deriv[signature]
	if verbose: print 'initial_guess', initial_guess
        deriv_val[signature], err_deriv_val[signature] = solve_relation_num(
	    drel[signature],
	    subsd,
	    deriv[signature],
	    initial_guess,
	    use_numexpr,
	    unitless,
	    bincompile,
	    verbose,
	    **kwargs
	    )
        subsd.update({deriv[signature]: deriv_val[signature]})

    return deriv_val, err_deriv_val


def test_solve_relation_for_derivatives():
    from sympy import symbols, Function, ln

    x, y = symbols('x, y')
    f = symbols('f', cls=Function)(x, y)
    relation = ln(x*y+f)-f
    subsd = {x:2,y:1}
    diff_wrt = {x:2,y:1}
    #diff_wrt_map = {x_arg: x, y_arg: y}
    df, df_err = solve_relation_for_derivatives(relation, {x:2,y:1},
						f, {((x,0),(y,0)):1.0}, diff_wrt,
						#diff_wrt_map,
						verbose=True,
						abstol=1e-8)
    assert abs(df[((x,1),(y,0))]-0.46594127238)       < 1e-8
    assert abs(df[((x,0),(y,1))]-0.9318825447699826)  < 1e-8
    assert abs(df[((x,1),(y,1))]+0.17057414955751987) < 1e-8
    assert abs(df[((x,2),(y,1))]-0.37120710736327495) < 1e-8

    x, y = symbols('x, y')
    #x_arg, y_arg = symbols('x_arg, y_arg')
    g = symbols('g', cls=Function)(x,y)#(x_arg, y_arg)
    relation = x**2+y**2-g # Circle
    subsd = {x: 1.0, y: 0.0}
    initial_guess_func_vals = {((x, 0), (y, 0)): 0.5}
    diff_wrt = {x: 1}
    dg, dg_err = solve_relation_for_derivatives(relation, subsd, g, initial_guess_func_vals, diff_wrt, verbose=True, abstol=1e-7)
    assert abs(dg[((x, 1),)]-2.0) < 1e-7


def get_cb_from_rel(rel, subsd, varied, use_numexpr=False, unitless=False, bincompile=False):
    """
    Turn a symbolic (sympy) equation into a python callback
    of one variable.

    Some rules of thumb:

    If the callback will be evaluated for large vectors set
        use_numexpr = True

    If the expression is very complicated and will be called pointwise set
        bincompile = True
    """
    from sympy import symbols
    from sympy.abc import x
    # If varied is general function dep. on (x,y) we cannot
    # substitute for x and y without error. Therefore we
    # mask `varied` as dummy and resubstitute it later
    dummy = symbols('_dummy_')
    subsrel = rel.subs({varied:dummy})
    subsrel = subsrel.subs(subsd)
    if bincompile:
	assert not use_numexpr
	assert unitless
	return autowrap(subsrel.subs({dummy:x}))
    else:
	if use_numexpr:
	    def f0(x):
		return ne.evaluate(str(subsrel), {str(dummy): x})
	else:
	    subsrel = subsrel.subs({dummy:varied})
	    def f0(x):
		return subsrel.subs({varied: x})
	return f0


def test_get_cb_from_rel():
    from sympy import symbols, Function, ln
    x,y = symbols('x,y')
    f = symbols('f',cls=Function)(x,y)
    relation = ln(x*y+f)-f
    subsd = {x:2,y:1}
    varied = f
    initial_guess = 1.0
    cb = get_cb_from_rel(relation, subsd, f)
    assert cb(1.14619322062) < 1e-12

    cb_numexpr = get_cb_from_rel(relation,subsd,f,use_numexpr=True)
    assert cb_numexpr(1.14619322062) < 1e-12

    cb_autowrap = get_cb_from_rel(relation,subsd,f,bincompile=True,unitless=True)
    assert cb_autowrap(1.14619322062) < 1e-12
