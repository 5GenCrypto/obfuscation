#!/usr/bin/env sage -python

from __future__ import print_function

from gradedencoding import GradedEncoding
from branchingprogram import (BranchingProgram, MATRIX_LENGTH)
import utils

from sage.all import flatten, Integer, load, MatrixSpace, ZZ

import collections, os, sys, time

MS = MatrixSpace(ZZ, MATRIX_LENGTH)

def ms2list(m):
    '''Convert an element in MS to a flat integer list'''
    m = [[int(e) for e in row] for row in m]
    return [int(e) for e in flatten(m)]

ObfLayer = collections.namedtuple('ObfLayer', ['inp', 'I', 'J'])

def load_obf(directory, inp, I, J):
    inp = load('%s/%s' % (directory, inp))
    I = load('%s/%s' % (directory, I))
    J = load('%s/%s' % (directory, J))
    return ObfLayer(int(inp), I, J)

def save_obf(layer, directory, idx):
    Integer(layer.inp).save('%s/%d.input' % (directory, idx))
    layer.I.save('%s/%d.I' % (directory, idx))
    layer.J.save('%s/%d.J' % (directory, idx))

class Obfuscator(object):
    def __init__(self, secparam, verbose=False):
        self.ge = GradedEncoding(verbose=verbose)
        self.secparam = secparam
        self.obfuscation = None
        self._verbose = verbose
        self.logger = utils.make_logger(self._verbose)
        self.logger('Obfuscation parameters:')
        self.logger('  Security Parameter: %d' % self.secparam)

    def save(self, directory):
        assert self.obfuscation is not None
        if not os.path.exists(directory):
            os.mkdir(directory)
        Integer(self.ge.x0).save('%s/x0' % directory)
        Integer(self.ge.pzt).save('%s/pzt' % directory)
        for idx, layer in enumerate(self.obfuscation):
            save_obf(layer, directory, idx)

    def load(self, directory):
        assert self.obfuscation is None
        x0 = load('%s/x0.sobj' % directory)
        pzt = load('%s/pzt.sobj' % directory)
        files = os.listdir(directory)
        inputs = sorted(filter(lambda s: 'input' in s, files))
        Is = sorted(filter(lambda s: 'I' in s, files))
        Js = sorted(filter(lambda s: 'J' in s, files))
        # XXX: will the order be preserved for >= 10 layers?
        self.obfuscation = [load_obf(directory, inp, I, J) for inp, I, J in
                            zip(inputs, Is, Js)]
        self.ge.load_system_params(self.secparam, len(self.obfuscation), x0,
                                   pzt)

    def _obfuscate_layer(self, layer):
        self.logger('Obfuscating layer...')
        start = time.time()
        m = ms2list(layer.I)
        m.extend(ms2list(layer.J))
        half = len(m) / 2
        es = self.ge.encode_list(m)
        I, J = MS(es[:half]), MS(es[half:])
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))
        return ObfLayer(layer.inp, I, J)

    def obfuscate(self, bp):
        self.ge.gen_system_params(self.secparam, len(bp))
        self.logger('Obfuscating...')
        start = time.time()
        self.obfuscation = [self._obfuscate_layer(layer) for layer in bp]
        end = time.time()
        self.logger('Obfuscation took: %f seconds' % (end - start))

    def evaluate(self, inp):
        assert self.obfuscation is not None
        comp = MS.identity_matrix()
        for m in self.obfuscation:
            comp = comp * (m.I if inp[m.inp] == '0' else m.J)
        return 1 if self.ge.is_zero(comp[0][0]) else 0
