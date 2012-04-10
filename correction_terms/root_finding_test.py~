#/usr/bin/env python
import unittest
from secant_method import *
from math import sin

# Define some test functions to use:

def parabola(x):
    # x2-x+0.244449
    # x = 0.5 pm sqrt(1e-6)
    return x*(x-1)+0.25-1e-6

parabola.roots = (0.5 - 1e-3, 0.5 + 1e-3)
parabola.start_guesses = (0, 1)


def poly3(x):
    return (x-0.9)*(x-1)*(x-1.1)
poly3.roots=(0.9,1.0,1.1)
poly3.start_guesses = (0.8,0.97,1.2)

# Test case

class TestSecantMethod(unittest.TestCase):

    def yabstol(self, func, x0, xtrue, dx0, tol, verbose):
        print "Testing yabstol"
        x = secant_method(func, x0, dx0, yabstol=tol, maxiter=30, verbose=verbose)
        print "f({0:17.10g}) = {1:17.10g}".format(x,func(x))
        print "analytical: {0:17.10g}, diff: {1:17.10g}".format(xtrue, x - xtrue)
        self.assertTrue(abs(func(x))<tol)

    def xabstol(self, func, x0, xtrue, dx0, tol, verbose):
        print "Testing xabstol"
        x = secant_method(func, x0, dx0, xabstol=tol, maxiter=30, verbose=verbose)
        print "f({0:20.13g}) = {1:17.10g}".format(x,func(x))
        print "analytical: {0:20.13g}, diff: {1:17.10g}".format(xtrue, x-xtrue)
        self.assertTrue(abs(x-xtrue)<tol)

    def abstol(self, func, x0, xtrue, dx0, tol, verbose):
        print "Testing abstol (xabstol and yabstol)"
        x = secant_method(func, x0, dx0, abstol=tol, maxiter=30, verbose=verbose)
        print "f({0:20.13g}) = {1:17.10g}".format(x,func(x))
        fmtstr="analytical: {0:"+self.precision+"}, diff: {1:"+ \
               self.precision+"}"
        print fmtstr.format(xtrue, x - xtrue)
        self.assertTrue(abs(x - xtrue)<tol and abs(func(x))<tol)
        

    def runTest(self):
        tol = 1e-6
        verbose = True
        self.precision = "20.13g"
        
        for func in [parabola, poly3]:
            for x0, xtrue in zip(func.start_guesses, func.roots):
                dx0 = (xtrue-x0)*1e-2
                self.yabstol(func, x0, xtrue, dx0, tol, verbose)
                self.xabstol(func, x0, xtrue, dx0, tol, verbose)
                self.abstol( func, x0, xtrue, dx0, tol, verbose)

if __name__ == '__main__':
    unittest.main()
