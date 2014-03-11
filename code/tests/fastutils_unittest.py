#!/usr/bin/env python2

from __future__ import print_function

import unittest

kappa = 1
alpha = 16
beta = 8
eta = 80
nu = 29
rho = 16
n = 80

class FastutilsTestCase(unittest.TestCase):
    def testme(self):
        try:
            import fastutils
            x0, pzt = fastutils.genparams(n, alpha, beta, eta, kappa)
            zero = fastutils.encode(long(0), rho)
            one = fastutils.encode(long(1), rho)
            assert not fastutils.is_zero(one, nu)
            assert fastutils.is_zero(zero, nu)
        except ImportError as e:
            self.fail(str(e))

def suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTest(loader.loadTestsFromTestCase(FastutilsTestCase))
    return suite
