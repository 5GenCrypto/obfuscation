#!/usr/bin/env python2

from __future__ import print_function

import numpy
import unittest

# kappa = 1
# alpha = 8
# beta = alpha
# eta = 58
# nu = 21
# rho = alpha
# n = eta

# kappa = 2
# alpha = 8
# beta = alpha
# eta = 76
# nu = 21
# rho = alpha
# n = eta

kappa = 3
alpha = 8
beta = alpha
eta = 94
nu = 21
rho = alpha
n = eta

g0 = long(181)

class FastutilsTestCase(unittest.TestCase):
    def testme(self):
        try:
            import fastutils

            nzs = 3

            # assert kappa == nzs

            x0, pzt = fastutils.genparams(n, alpha, beta, eta, kappa, nzs, g0)

            # p = 143L
            s = [55L, 103L, 50L, 93L, 151L]
            t = [59L, 180L, 57L, 168L, 34L]
            p = long(numpy.dot(s, t) % g0)
            assert p == 143L

            M = numpy.matrix([[1L, 0L, 0L, 0L, 0L],
                              [0L, 1L, 0L, 0L, 0L],
                              [0L, 0L, 1L, 0L, 0L],
                              [0L, 0L, 0L, 1L, 0L],
                              [0L, 0L, 0L, 0L, 1L]])

            l = M.flatten().tolist()[0]
            Menc = fastutils.encode_vector(l, rho, 0)
            Menc = numpy.array(Menc).reshape((5, 5))

            print(type(Menc))

            penc = fastutils.encode_scalar(p, rho, 1, 2)
            senc = fastutils.encode_vector(s, rho, 1)
            tenc = fastutils.encode_vector(t, rho, 2)
            one = fastutils.encode_scalar(1L, rho, 0, -1)

            a = numpy.dot(numpy.dot(senc, Menc), tenc) % x0
            b = (penc * fastutils.encode_scalar(1L, rho, 0, -1)) % x0

            test = a - b

            print()
            print('test length    = %d' % long(test).bit_length())
            assert fastutils.is_zero(long(test), nu)
        except ImportError as e:
            self.fail(str(e))

def suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTest(loader.loadTestsFromTestCase(FastutilsTestCase))
    return suite
