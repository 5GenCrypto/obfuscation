#!/usr/bin/env python2

from __future__ import print_function
import unittest
import doctest

alpha = 16
beta = 8
eta = 360
kappa = 4
rho = 16
n = 360

ms = [0] * n
ms[0] = 1

class DeviceTest(unittest.TestCase):
    def runTest(self):
        try:
            import fastutils
            x0, pzt = fastutils.genparams(n, alpha, beta, eta, kappa)
        except ImportError as e:
            self.fail(str(e))

