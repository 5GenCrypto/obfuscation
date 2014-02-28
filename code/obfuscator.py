#!/usr/bin/env sage -python

from __future__ import print_function

from gradedencoding import GradedEncoding
from branchingprogram import (BranchingProgram, MATRIX_LENGTH)
import utils

from sage.all import *

import functools, time, sys

MS = MatrixSpace(ZZ, MATRIX_LENGTH)

class ObfLayer(object):
    def __init__(self, inp, I, J):
        self.inp = inp
        self.I = I
        self.J = J
    def __repr__(self):
        return "%d\n%s\n%s" % (self.inp, self.I, self.J)

class Obfuscator(object):
    def __init__(self, secparam, bp, verbose=False, parallel=False, ncpus=1):
        self.ge = GradedEncoding(secparam, len(bp), verbose=verbose,
                                 parallel=parallel, ncpus=ncpus)
        self.obfuscation = None
        self.bp = bp
        self._verbose = verbose
        self._parallel = parallel
        self.logger = functools.partial(utils.logger, verbose=self._verbose)

    def load(self, fname):
        raise NotImplemented()

    def _obfuscate_matrix(self, m):
        m = [[int(e) for e in row] for row in m]
        m = [int(e) for e in flatten(m)]
        self.logger('Obfuscating matrix %s' % m)
        start = time.time()
        if self._parallel:
            m = self.ge.encode_list(m)
        else:
            m = [self.ge.encode(e) for e in m]
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))
        return MS(m)

    def _obfuscate_layer(self, layer):
        I = self._obfuscate_matrix(layer.I)
        J = self._obfuscate_matrix(layer.J)
        return ObfLayer(layer.inp, I, J)

    def obfuscate(self):
        self.obfuscation = [self._obfuscate_layer(layer) for layer in self.bp]

    def save(self, fname):
        assert self.obfuscation is not None
        raise NotImplemented()

    def _mult_matrices(self, A, B):
        rows = []
        for row in A:
            elems = []
            for col in B.transpose():
                m = [self.ge.mult((row[i], col[i])) for i in xrange(MATRIX_LENGTH)]
                a = self.ge.add(m)
                elems.append(a)
            rows.append(elems)
        return MS(rows)

    def evaluate(self, inp):
        assert self.obfuscation is not None
        comp = MS.identity_matrix()
        for m in self.obfuscation:
            comp = self._mult_matrices(comp, m.I if inp[m.inp] == '0' else m.J)
        return 1 if self.ge.is_zero(comp[0][0]) else 0
