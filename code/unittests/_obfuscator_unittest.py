
from __future__ import print_function

import numpy as np
import unittest

from indobf.branchingprogram import BranchingProgram, ParseException
from indobf.obfuscator import Obfuscator
import indobf.test as test
import indobf._obfuscator as _obf

params = test.TestParams(obfuscate=True)
verbose = False

class EncodeTestCase(unittest.TestCase):

    # def test_onelevel(self):
    #     n = 174
    #     alpha = 8
    #     beta = alpha
    #     eta = 58
    #     nu = 21
    #     rho = alpha
    #     nzs = 1
    #     gs = [long(181)]

    #     zero = np.array([[1L, 0L, 0L, 0L, 0L, 0L],
    #                      [0L, 1L, 0L, 0L, 0L, 0L],
    #                      [0L, 0L, 1L, 0L, 0L, 0L],
    #                      [0L, 0L, 0L, 1L, 0L, 0L],
    #                      [0L, 0L, 0L, 0L, 1L, 0L],
    #                      [0L, 0L, 0L, 0L, 0L, 1L]], dtype=object)
    #     one = np.array([[0L, 0L, 0L, 1L, 0L, 0L],
    #                     [0L, 0L, 1L, 0L, 0L, 0L],
    #                     [0L, 1L, 0L, 0L, 0L, 0L],
    #                     [1L, 0L, 0L, 0L, 0L, 0L],
    #                     [0L, 0L, 0L, 0L, 1L, 0L],
    #                     [0L, 0L, 0L, 0L, 0L, 1L]], dtype=object)

    #     _obf.setup(n, alpha, beta, eta, nu, rho, nzs, gs, 'testdir')
    #     _obf.encode_level(0, 0, zero, one, [0], [0])
    #     assert _obf.evaluate('testdir', '0', 1) == 0
    #     assert _obf.evaluate('testdir', '1', 1) == 1

    def test_id_circuit_8(self):
        assert test.test_circuit('circuits/id.circ', 8, verbose, params)

    # def test_not_circuit_8(self):
    #     assert test.test_circuit('circuits/not.circ', 8, verbose, params)

    # def test_xor_circuit_8(self):
    #     assert test.test_circuit('circuits/xor.circ', 8, verbose, params)

    # def test_and_circuit_8(self):
    #     assert test.test_circuit('circuits/and.circ', 8, verbose, params)

    def test_id_circuit_16(self):
        assert test.test_circuit('circuits/id.circ', 16, verbose, params)

    # def test_id_circuit_32(self):
    #     assert test.test_circuit('circuits/id.circ', 32, verbose, params)
