#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This work is open source and is released under the
# 2-clause BSD license (see LICENSE.txt for further information)
# Copyright (c) 2011, 2012, Björn Ingvar Dahlgren

"""
Currently not used in the calculation of correction terms.
"""

import numpy as np

def get_derivative_from_finite_difference(order,x_list,y_list,x0):
    """
    Algorithm from:
    Generation of Finite Difference Formulas on Arbitrarily Spaced Grids, Bengt Fornberg
    Mathematics of compuation, 51, 184, 1988, 699-706
    """
    if len(x_list) != len(y_list):
        raise NameError("x_list and y_list not equal in length.")

    M=order; N=len(x_list)-1
    delta = np.zeros((M+1,N+1,N+1))
    delta[0,0,0]=1.0
    c1 = 1.0
    for n in range(1,N+1):
        c2 = 1.0
        for nu in range(0,n):
            c3 = x_list[n]-x_list[nu]
            c2 = c2 * c3
            if n <= M: delta[n,n-1,nu]=0
            for m in range(0,min(n,M)+1):
                delta[m,n,nu] = ( (x_list[n]-x0)*delta[m,n-1,nu] - m*delta[m-1,n-1,nu] )/c3
        for m in range(0,min(n,M)+1):
            delta[m,n,n] = c1/c2*( m*delta[m-1,n-1,n-1] - (x_list[n-1]-x0)*delta[m,n-1,n-1] )
        c1 = c2
    derivative = 0
    for nu in range(0,N+1):
        derivative += delta[M,N,nu]*y_list[nu]
    return derivative


def test():
    n = 20
    l = 0.25
    ii = 0
    repl_lst = []
    while ii < 10:
        ii += 1
        x = np.cumsum(np.random.rand(n)*l)-0.25*n*l
        #print x
        func = lambda x: np.exp(x)
        lst = []
        for order in range(4):
            lst.append(get_derivative_from_finite_difference(order,x,func(x),0))
        if np.any(np.array(lst)>2):
            print x
        repl_lst.append(lst)
    print np.array(repl_lst)

def test_err():
    n = 20
    l = 0.25
    ii = 0
    repl_lst = []
    while ii < 3:
        ii += 1
        x = np.cumsum(np.random.rand(n)*l)
        x = np.concatenate((x[::-1]*-1,x))
        #print x
        func = lambda x: np.exp(x)
        lst = []
        derivs = np.zeros(n)
        for order in range(4):
            print 'Order',order
            for m in range(1,n+1):
                sub_x=x[n-m:n+m]
                derivs[m-1] = get_derivative_from_finite_difference(order,sub_x,func(sub_x),0)
            print '{0: <5s}\t{1: <21s}\t{2: >21s}\t{3: >21s}\t{4: >21s}'.format('m','val','diff','analytical error','diff/analytical')
            for m in range(1,n):
                print '{0: <5d}\t{1:20.18f}\t{2: >21.18f}\t{3: >21.18f}\t{4: >21.18f}'.format((m+1)*2,derivs[m],derivs[m]-derivs[m-1],derivs[m]-1,(derivs[m]-derivs[m-1])/(derivs[m]-1))
            lst.append((derivs[-1],abs(derivs[-1]-derivs[-2])))
        repl_lst.append(lst)
    print np.array(repl_lst)
