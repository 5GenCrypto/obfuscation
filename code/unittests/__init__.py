#!/usr/bin/env python2

import unittest

def suite():
    suite = unittest.TestSuite()
    suite.addTests(fastutils_unittest.suite())
    return suite

if __name__ == '__main__':
    unittest.TextTestRunner(verbosity=2).run(suite())
