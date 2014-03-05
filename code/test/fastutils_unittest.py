#!/usr/bin/env python2

from __future__ import print_function
import unittest
import doctest

class DeviceTest(unittest.TestCase):
    def runTest(self):
        try:
            import fastutils
            primes = fastutils.genprimes(100, 100)
            print(primes)
        except ImportError as e:
            self.fail(str(e))

