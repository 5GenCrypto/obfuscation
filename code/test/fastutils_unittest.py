#!/usr/bin/env python2

from __future__ import print_function
import unittest
import doctest

alpha = 16
beta = 8
eta = 360
kappa = 4
n = 360

class DeviceTest(unittest.TestCase):
    def runTest(self):
        try:
            import fastutils
            x0, ps, gs, z, zinv, pzt \
              = fastutils.genparams(n, alpha, beta, eta, kappa)
            # print(x0)
            # print(ps)
            # print(gs)
            # print(z)
            # print(zinv)
            # print(pzt)
            # fastutils.init(num, ps, gs, 0, 0)
            # primes = fastutils.genprimes(100, 100)
            # x0, primes = fastutils.genprimes(100, 100)
            # print(primes)
        except ImportError as e:
            self.fail(str(e))

