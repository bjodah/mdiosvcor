#/usr/bin/env python
from __future__ import division
import unittest
from root_finding import *
from math import sin, cos

# Define some test functions to use:

def parabola(x):
    # x2-x+0.244449
    # x = 0.5 pm sqrt(1e-8)
    return x*(x-1)+0.25-1e-8

def parabola_der(x):
    return x-1+x

parabola.roots         = (0.5 - 1e-4, 0.5 + 1e-4)
parabola.start_guesses = (0, 1)
parabola.dfdx          = parabola_der
parabola.name          = "f(x) = x*(x-1)+0.25-1e-8"

def poly3(x):
    return (x-0.9)*(x-1)*(x-1.1)

def poly3_der(x):
    return (x-1)*(x-1.1)+(x-0.9)*(x-1)+(x-0.9)*(x-1.1)

poly3.roots         = (0.9,1.0,1.1)
poly3.start_guesses = (0.8,0.97,1.2)
poly3.dfdx          = poly3_der
poly3.name          = "f(x) = (x-0.9)*(x-1)*(x-1.1)"

# Test case

def ex3_3(x):
    """ From P. Pohl p. 77"""
    return x - 4*sin(2*x)-3-3/80

def ex3_3_der(x):
    return 1 - 8*cos(2*x)

def ex3_3_fix_point(x):
    """ Ex 3.5 NR modified P. Pohl p. 87"""
    return x-(x-4*sin(2*x)-3-3/80)/9

ex3_3.dfdx = ex3_3_der
ex3_3.fix_point = ex3_3_fix_point

class TestSecantMethod(unittest.TestCase):
    def runTest(self):
        x, y = find_root(parabola, parabola.start_guesses[0], 'secant',
                         parabola.dfdx, verbose = True)
        

class TestTolerances(unittest.TestCase):

    def yabstol(self, func, x0, preferred_method, dfdx, xtrue, dx0, tol, verbose):
        print "Testing yabstol"
        x, y = find_root(func, x0, preferred_method, dfdx, dx0, yabstol=tol, maxiter=30, verbose=verbose)
        print "f({0:20.13g}) = {1: <17.10g}".format(x,func(x))
        fmtstr="analytical: {0: <"+self.precision+"}, diff: {1: <"+ \
               self.precision+"}"
        print fmtstr.format(xtrue, x - xtrue)
        self.assertTrue(abs(func(x))<tol)

    def xabstol(self, func, x0, preferred_method, dfdx, xtrue, dx0, tol, verbose):
        print "Testing xabstol"
        x, y = find_root(func, x0, preferred_method, dfdx, dx0, xabstol=tol, maxiter=30, verbose=verbose)
        print "f({0:20.13g}) = {1: <17.10g}".format(x,func(x))
        fmtstr="analytical: {0: <"+self.precision+"}, diff: {1: <"+ \
               self.precision+"}"
        print fmtstr.format(xtrue, x - xtrue)
        self.assertTrue(abs(x-xtrue)<tol)

    def abstol(self, func, x0, preferred_method, dfdx, xtrue, dx0, tol, verbose):
        print "Testing abstol (xabstol and yabstol)"
        x, y = find_root(func, x0, preferred_method, dfdx, dx0, abstol=tol, maxiter=30, verbose=verbose)
        print "f({0:20.13g}) = {1: <17.10g}".format(x,func(x))
        fmtstr="analytical: {0: <"+self.precision+"}, diff: {1: <"+ \
               self.precision+"}"
        print fmtstr.format(xtrue, x - xtrue)
        self.assertTrue(abs(x - xtrue)<tol and abs(func(x))<tol)
        

    def runTest(self):
        tol = 1e-9
        verbose = False
        self.precision = "20.13g"
        for preferred_method in ["newton","secant"]:
            for func in [parabola, poly3]:
                print func.name
                dfdx = func.dfdx
                for x0, xtrue in zip(func.start_guesses, func.roots):
                    dx0 = (xtrue-x0)*1e-2
                    self.yabstol(func, x0, preferred_method, dfdx, xtrue, dx0, tol, verbose)
                    self.xabstol(func, x0, preferred_method, dfdx, xtrue, dx0, tol, verbose)
                    self.abstol( func, x0, preferred_method, dfdx, xtrue, dx0, tol, verbose)

if __name__ == '__main__':
    unittest.main()
