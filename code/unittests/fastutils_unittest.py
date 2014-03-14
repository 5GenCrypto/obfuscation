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

nzs = 1

g0 = long(129)

class FastutilsTestCase(unittest.TestCase):
    def testme(self):
        try:
            import fastutils
            x0, pzt = fastutils.genparams(n, alpha, beta, eta, kappa, nzs, g0)
            a = fastutils.encode_scalar(long(1), rho, 0, -1)
            b = fastutils.encode_scalar(long(1), rho, 0, -1)
            # c = fastutils.encode_scalar(long(1), rho, 0, 1)
            # test = a * b - c
            test = abs(a - b)
            print(test)
            assert fastutils.is_zero(test, nu)
            test = abs(b - a)
            print(test)
            assert fastutils.is_zero(test, nu)
        except ImportError as e:
            self.fail(str(e))

def suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTest(loader.loadTestsFromTestCase(FastutilsTestCase))
    return suite
