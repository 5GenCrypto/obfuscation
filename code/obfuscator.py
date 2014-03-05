#!/usr/bin/env sage -python

from __future__ import print_function

from gradedencoding import GradedEncoding
from branchingprogram import (BranchingProgram, MATRIX_LENGTH)
import utils

from sage.all import *

import os, sys, time

MS = MatrixSpace(ZZ, MATRIX_LENGTH)

def ms2list(m):
    '''Convert an element in MS to a flat integer list'''
    m = [[int(e) for e in row] for row in m]
    return [int(e) for e in flatten(m)]

class ObfLayer(object):
    def __init__(self, inp, I, J):
        self.inp = inp
        self.I = I
        self.J = J
    @classmethod
    def load(cls, directory, inp, I, J):
        inp = load('%s/%s' % (directory, inp))
        I = load('%s/%s' % (directory, I))
        J = load('%s/%s' % (directory, J))
        return cls(int(inp), I, J)
    def save(self, directory, idx):
        Integer(self.inp).save('%s/%d.input' % (directory, idx))
        self.I.save('%s/%d.I' % (directory, idx))
        self.J.save('%s/%d.J' % (directory, idx))
    def __repr__(self):
        return "%d\n%s\n%s" % (self.inp, self.I, self.J)

class Obfuscator(object):
    def __init__(self, secparam, verbose=False, parallel=False, ncpus=1,
                 use_c=False):
        self.ge = GradedEncoding(verbose=verbose, parallel=parallel,
                                 ncpus=ncpus, use_c=use_c)
        self.secparam = secparam
        self.obfuscation = None
        self._verbose = verbose
        self._parallel = parallel
        self._use_c = use_c
        self.logger = utils.make_logger(self._verbose)
        self.logger('Obfuscation parameters:')
        self.logger('  Security Parameter: %d' % self.secparam)
        self.logger('  Parallel: %s' % self._parallel)
        self.logger('  Verbose: %s' % self._verbose)
        self.logger('  Using C: %s' % self._use_c)

    def load(self, directory):
        assert self.obfuscation is None
        x0 = load('%s/x0.sobj' % directory)
        pzt = load('%s/pzt.sobj' % directory)
        # REFACTOR: this is a mess
        files = os.listdir(directory)
        files.remove('x0.sobj')
        files.remove('pzt.sobj')
        inputs = filter(lambda s: 'input' in s, files)
        inputs.sort()
        Is = filter(lambda s: 'I' in s, files)
        Is.sort()
        Js = filter(lambda s: 'J' in s, files)
        Js.sort()
        # XXX: will the order be preserved for >= 10 layers?
        self.obfuscation = [ObfLayer.load(directory, inp, I, J) for inp, I, J in
                            zip(inputs, Is, Js)]
        self.ge.load_system_params(self.secparam, len(self.obfuscation), x0,
                                   pzt)

    def _obfuscate_matrix(self, m):
        m = ms2list(m)
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

    def obfuscate(self, bp):
        self.ge.gen_system_params(self.secparam, len(bp))
        self.obfuscation = [self._obfuscate_layer(layer) for layer in bp]

    def save(self, directory):
        assert self.obfuscation is not None
        if not os.path.exists(directory):
            os.mkdir(directory)
        self.ge.x0.save('%s/x0' % directory)
        self.ge.pzt.save('%s/pzt' % directory)
        for idx, layer in enumerate(self.obfuscation):
            layer.save(directory, idx)

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
