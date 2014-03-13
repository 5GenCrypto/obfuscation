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

ObfLayer = collections.namedtuple('ObfLayer', ['inp', 'zero', 'one'])

def load_obf(directory, inp, zero, one):
    inp = load('%s/%s' % (directory, inp))
    zero = load('%s/%s' % (directory, zero))
    one = load('%s/%s' % (directory, one))
    return ObfLayer(int(inp), zero, one)

def save_obf(layer, directory, idx):
    Integer(layer.inp).save('%s/%d.input' % (directory, idx))
    layer.zero.save('%s/%d.zero' % (directory, idx))
    layer.one.save('%s/%d.one' % (directory, idx))

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
        zeros = sorted(filter(lambda s: 'zero' in s, files))
        ones = sorted(filter(lambda s: 'one' in s, files))
        # XXX: will the order be preserved for >= 10 layers?
        self.obfuscation = [load_obf(directory, inp, zero, one) for inp, zero,
                            one in zip(inputs, zeros, ones)]
        self.ge.load_system_params(self.secparam, len(self.obfuscation), x0,
                                   pzt)

    def _obfuscate_layer(self, layer, idx):
        self.logger('Obfuscating layer...')
        start = time.time()
        m = ms2list(layer.zero)
        m.extend(ms2list(layer.one))
        half = len(m) / 2
        es = self.ge.encode_layer(m, idx)
        zero, one = MS(es[:half]), MS(es[half:])
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))
        return ObfLayer(layer.inp, zero, one)

    def obfuscate(self, bp, g0=None):
        self.ge.gen_system_params(self.secparam, len(bp), g0=g0)
        self.logger('Obfuscating...')
        start = time.time()
        self.obfuscation = [self._obfuscate_layer(layer, idx) for idx, layer in
                            enumerate(bp)]
        end = time.time()
        self.logger('Obfuscation took: %f seconds' % (end - start))

    def evaluate(self, inp):
        assert self.obfuscation is not None
        comp = MS.identity_matrix()
        for m in self.obfuscation:
            comp = comp * (m.zero if inp[m.inp] == '0' else m.one)
        if self.ge.is_zero(comp[0][1]) and self.ge.is_zero(comp[1][0]):
            return 0
        else:
            return 1
