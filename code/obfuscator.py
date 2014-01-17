#!/usr/bin/env sage -python

from __future__ import print_function

from gradedencoding import GradedEncoding
from branchingprogram import (BranchingProgram, MATRIX_LENGTH)
from sage.all import *
import sys

MS = MatrixSpace(ZZ, MATRIX_LENGTH)

class ObfLayer(object):
    def __init__(self, inp, I, J):
        self.inp = inp
        self.I = I
        self.J = J
    def __repr__(self):
        return "%d\n%s\n%s" % (self.inp, self.I, self.J)

class Obfuscator(object):
    def __init__(self, secparam, verbose=False):
        self.ge = GradedEncoding(secparam, len(bp), verbose=verbose)
        self.obfuscation = None
        self._verbose = verbose

    def load(self, fname):
        raise NotImplemented()

    def _obfuscate_matrix(self, m):
        rows = []
        for i, row in enumerate(m):
            elems = []
            for j, elem in enumerate(row):
                r = self.ge.encode(int(elem))
                elems.append(r)
            rows.append(elems)
        return MS(rows)

    def _obfuscate_layer(self, layer):
        I = self._obfuscate_matrix(layer.I)
        J = self._obfuscate_matrix(layer.J)
        return ObfLayer(layer.inp, I, J)

    def obfuscate(self, bp):
        self.obfuscation = [self._obfuscate_layer(layer) for layer in bp]

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
        assert(self.obfuscation is not None)
        comp = MS.identity_matrix()
        for m in self.obfuscation:
            comp = self._mult_matrices(comp, m.I if inp[m.inp] == '0' else m.J)
        return 1 if self.ge.is_zero(comp[0][0]) else 0

if __name__ == '__main__':
    fname = sys.argv[1]
    inp = sys.argv[2]
    print('Converting circuit -> bp...')
    bp = BranchingProgram(fname, type='circuit')
    # bp.obliviate()
    bp.randomize()
    print(bp)
    print('Obfuscating BP of length %d...' % len(bp))
    obf = Obfuscator(8, verbose=True)
    obf.obfuscate(bp)
    print('Evaluating on input %s...' % inp)
    r = obf.evaluate(inp)
    print('Output = %d' % r)
