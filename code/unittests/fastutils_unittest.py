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

gs = [long(181)]

class FastutilsTestCase(unittest.TestCase):
    def testme(self):
        try:
            import indobf.fastutils as fastutils

            nzs = 3

            # assert kappa == nzs

            x0, pzt = fastutils.genparams(n, alpha, beta, eta, kappa, rho, nzs, gs)

            # a1 = fastutils.encode_scalar([173L, 173L], [0, 1], [0, 1]);
            # b1 = fastutils.encode_scalar([173L, 173L], [0, 1], [0, 1]);

            # assert fastutils.is_zero(long(a1 - b1), nu)

            p1 = 143L
            s1 = [55L, 103L, 50L, 93L, 151L]
            t1 = [59L, 180L, 57L, 168L, 34L]

            M = numpy.matrix([[1L, 0L, 0L, 0L, 0L],
                              [0L, 1L, 0L, 0L, 0L],
                              [0L, 0L, 1L, 0L, 0L],
                              [0L, 0L, 0L, 1L, 0L],
                              [0L, 0L, 0L, 0L, 1L]])

            l = M.flatten().tolist()[0]
            M1enc = fastutils.encode_vector(l, 0, [0])
            M1enc = numpy.array(M1enc).reshape((5, 5))

            p1enc = fastutils.encode_scalar(p1, 0, [1, 2])
            s1enc = fastutils.encode_vector(s1, 0, [1])
            t1enc = fastutils.encode_vector(t1, 0, [2])

            a = numpy.dot(numpy.dot(s1enc, M1enc), t1enc) % x0
            b = (p1enc * fastutils.encode_scalar(1L, 0, [0])) % x0

            assert fastutils.is_zero(long(a - b), nu)

        except ImportError as e:
            self.fail(str(e))

def suite():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTest(loader.loadTestsFromTestCase(FastutilsTestCase))
    return suite
