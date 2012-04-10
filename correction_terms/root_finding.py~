#!/usr/bin/python
# -*- coding: utf-8 -*-

from functools import reduce
from operator import and_

def secant_method(f,
                  x0,
                  dx0=None,
                  unit=None,
                  abstol=None,
                  yabstol=1e-8,
                  xabstol=None,
                  maxiter=25,
                  return_intermediate_results=False,
                  verbose=False):
    def cb(x):
        if unit:
            return f(x*unit)
        else:
            return f(x)

    satisfication_functions = []

    def satisfied():
        res = [f() for f in satisfication_functions]
        return reduce(and_,res)

    def satisfied_yabstol():
        return abs(y)<yabstol

    def satisfied_xabstol():
        return abs(dx)<xabstol

    if abstol:
        yabstol = abstol
        xabstol = abstol
        
    if yabstol: satisfication_functions.append(satisfied_yabstol)
    if xabstol: satisfication_functions.append(satisfied_xabstol)

    i = 0
    if return_intermediate_results:
        iterm_res = []

    # First step in secant method
    y0 = cb(x0)
    if return_intermediate_results:
        iterm_res.append((x0,y0))

    if verbose:
        fmtstr="{0: >3}{1: >14}{2: >14}{3: >14}{4: >14}{5: >14}{6: >15}"
        print fmtstr.format("i","x","dx","y","dy","dxi/dx(i-1)","dxi2/sqrt(dx(i-1))")
        print fmtstr.format(i,x0,"-",y0,"-","-","-")

    if dx0:
        dx = dx0
    else:
        dx = x0*1e-2 + 1e-6 # arbitrary first step
    x  = x0+dx
    y  = cb(x)
    if return_intermediate_results:
        res.append((x,y))

    dy = y-y0
    i  += 1
    if verbose:
        print fmtstr.format(i,x,dx,y,dy,"-","-")
    while not satisfied() and i <= maxiter:
        if verbose: dx0 = dx
        dx = -y*dx/dy
        x0 = x
        x += dx
        y0 = y
        y  = cb(x)
        if return_intermediate_results:
            interm_res.append((x,y))

        dy = y-y0
        if verbose:
            print fmtstr.format(i,x,dx,y,dy,dx0/dx,dx0/abs(dx)**0.5)
        i += 1
        
        
    if i > maxiter:
        raise ValueError
    else:
        if return_intermediate_results:
            return interm_res
        else:
            return x

