#!/usr/bin/python
# -*- coding: utf-8 -*-

# Code for implementing common
# root finding algorithms, implemented
# with sympy.physics.units compability
# in mind. (I long for the day sympy
# interacts well with the PyPI package:
# "quantities")

# Author: Björn Dahlgren
# Implemented for use in research project in the IGC group at ETH

# To the extent possible under law, Bjoern Dahlgren has waived all
# copyright and related or neighboring rights to this work.


from functools import reduce
from operator import and_

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
                 fmtstr="{0: >3}{1: >25}{2: >25}"+\
                        "{3: >25}{4: >25}{5: >25}"):
    """
    For studying the convergence it is helpful to
    see how the values changes between steps.
    Call find_root with "verbose=True" in order
    to enable printing of intermediate values.
    """

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
        print fmtstr.format(i,x/xunit,"-",y/yunit,"-","-")
    elif i ==1:
        print fmtstr.format(i,x/xunit,dx/xunit,y/yunit,dy/yunit,'-')
    else:
        print fmtstr.format(i,x/xunit,dx/xunit,y/yunit,dy/yunit,dx/dx_old**2*xunit)

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
              maxiter=25,
              return_intermediate_results=False,
              verbose=False):
    """
    If no preferred method is specified:
        If dfdx is specified, assume that Newtons method is wanted.
        Else assume secant method.

    method choice is stored as a tuple:
    (callback, args)

    TODO: let unit be automatically taken from x0 and func
    """

    if preferred_method == "newton" or \
        (not preferred_method and dfdx):
        if verbose: print "Using Newton's method."
        method = (newton_generator, (func, x0, dfdx))
    elif preferred_method == "secant" or \
        (not preferred_method and not dfdx):
        if verbose: print "Using the secant method."
        if not dx0:
            # If no dx0 specified take an arbitrary step forward
            dx0 = x0*1e-2 + 1e-6*xunit
            if verbose: print "Assuming dx0 = {0}".format(dx0)
        method = (secant_generator, (func, x0, dx0))
    else:
        NotImplemented


    # Setup tolerance requierments
    satisfication_functions = []

    def satisfied():
        satisfication = [sf() for sf in satisfication_functions]
        return reduce(and_,satisfication)

    def satisfied_yabstol():
        return abs(y/yunit)<yabstol

    def satisfied_xabstol():
        return abs(dx/xunit)<xabstol

    if abstol:
        if xabstol or yabstol:
            raise ValueError("Both abstol and xabstol and/or yabstol specified")
        yabstol = abstol
        xabstol = abstol

    if not yabstol and not xabstol:
        yabstol = 1e-10
        if verbose: print "Since no tolerances were specified," + \
           " abstol assumed = {0:12.5g}".format(yabstol)

    if yabstol: satisfication_functions.append(satisfied_yabstol)
    if xabstol: satisfication_functions.append(satisfied_xabstol)

    # Handle itermediate results if requested
    if return_intermediate_results:
        interm_res = []

    def append_res(x,y):
	iterm_res.append((x,y))


    # Initialize the generator of chosen method
    method_gen = method[0](*method[1])
    i = 0
    x0, y0 = method_gen.next()
    dx, dx_old, dy = None, None, None
    if verbose: print_info(i,x0,dx,y0,dy,dx_old,xunit,yunit)
    for x,y in method_gen:
        i     += 1
        dx_old = dx
        dx, dy = x - x0, y - y0
        if verbose: print_info(i,x,dx,y,dy,dx_old,xunit,yunit)
        if return_intermediate_results: append_res(x,y)
        if satisfied(): break
        if i == maxiter:
            raise ValueError("Maximum number of iterations reached")
        x0, y0 = x, y

    if return_intermediate_results:
        return interm_res
    else:
	return x,y

