#!/usr/bin/python
# -*- coding: utf-8 -*-

from functools import reduce
from operator import and_

def secant_generator(f, x0, dx0):
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

def print_info(i,x,dx,y,dy,dx_old,
                 fmtstr="{0: >3}{1: >20}{2: >20}"+\
                        "{3: >20}{4: >20}{5: >20}"):
    if i==0:
        print fmtstr.format("i","x","dx","y","dy","dxi/dx(i-1)**2")
        print fmtstr.format(i,x,"-",y,"-","-")
    elif i ==1:
        print fmtstr.format(i,x,dx,y,dy,'-')
    else:
        print fmtstr.format(i,x,dx,y,dy,dx/dx_old**2)

def find_root(f,
              x0,
              preferred_method=None,
              dfdx=None,
              dx0=None,
              unit=None,
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
    """

    if preferred_method == "newton" or \
        (not preferred_method and dfdx):
        if verbose: print "Using Newton's method."
        method = (newton_generator, (f, x0, dfdx))
    elif preferred_method == "secant" or \
        (not preferred_method and not dfdx):
        if verbose: print "Using the secant method."
        if not dx0:
            # If no dx0 specified take an arbitrary step forward
            dx0 = x0*1e-2 + 1e-6
            if verbose: print "Assuming dx0 = {0:12.5g}".format(dx0)
        method = (secant_generator, (f, x0, dx0))
    else:
        NotImplemented

    # Handle units if x has got a unit
    if unit:
        try:
            float(x0/unit)
        except:
            ret_unit = None
        else:
            x0 = float(x0/unit)
            ret_unit = unit
        
    # Setup tolerance requierments
    satisfication_functions = []

    def satisfied():
        satisfication = [f() for f in satisfication_functions]
        return reduce(and_,satisfication)

    def satisfied_yabstol():
        return abs(y)<yabstol

    def satisfied_xabstol():
        return abs(dx)<xabstol

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
        if unit:
            interm_res.append((x*unit,f(x*unit)))
        else:
            interm_res.append((x,y))
        

    # Initialize the generator of chosen method
    method_gen = method[0](*method[1])
    i = 0
    x0, y0 = method_gen.next()
    dx, dx_old, dy = None, None, None
    if verbose: print_info(i,x0,dx,y0,dy,dx_old)
    for x,y in method_gen:
        i     += 1
        dx_old = dx
        dx, dy = x - x0, y - y0
        if verbose: print_info(i,x,dx,y,dy,dx_old)
        if return_intermediate_results: append_res(x,y)
        if satisfied(): break
        if i == maxiter:
            raise ValueError("Maximum number of iterations reached")
        x0, y0 = x, y

    if return_intermediate_results:
        return interm_res
    else:
        return x,y
    
