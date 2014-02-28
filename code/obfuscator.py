#!/usr/bin/env sage -python

from __future__ import print_function

from gradedencoding import GradedEncoding
from branchingprogram import (BranchingProgram, MATRIX_LENGTH)
import utils

from sage.all import *

import functools, json, time, sys

MS = MatrixSpace(ZZ, MATRIX_LENGTH)

def ms2list(m):
    m = [[int(e) for e in row] for row in m]
    return [int(e) for e in flatten(m)]

class ObfLayer(object):
    def __init__(self, inp, I, J):
        self.inp = inp
        self.I = I
        self.J = J
    def __repr__(self):
        return "%d\n%s\n%s" % (self.inp, self.I, self.J)

class ObfEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, ObfLayer):
            return [obj.inp, ms2list(obj.I), ms2list(obj.J)]
        return json.JSONEncoder.default(self, obj)

class Obfuscator(object):
    def __init__(self, secparam, verbose=False, parallel=False, ncpus=1):
        self.ge = GradedEncoding(verbose=verbose, parallel=parallel,
                                 ncpus=ncpus)
        self.secparam = secparam
        self.obfuscation = None
        self.bp = None
        self._verbose = verbose
        self._parallel = parallel
        self.logger = functools.partial(utils.logger, verbose=self._verbose)

    def load_bp(self, bp):
        self.ge.gen_system_params(self.secparam, len(bp))
        self.bp = bp

    def load_obf(self, directory):
        assert self.obfuscation is None
        assert os.path.exists(directory)
        self.obfuscation = []
        x0 = load('%s/x0.sobj' % directory)
        pzt = load('%s/pzt.sobj' % directory)
        self.ge.load_system_params(self.secparam, len(self.obfuscation), x0,
                                   pzt)
        files = os.listdir(directory)
        files.remove('x0.sobj')
        files.remove('pzt.sobj')
        inputs = filter(lambda s: 'input' in s, files)
        Is = filter(lambda s: 'I' in s, files)
        Js = filter(lambda s: 'J' in s, files)
        self.obfuscation = []
        # XXX: will the order be preserved for multiple layers?
        for inpfile, Ifile, Jfile in zip(inputs, Is, Js):
            inp = load('%s/%s' % (directory, inpfile))
            I = load('%s/%s' % (directory, Ifile))
            J = load('%s/%s' % (directory, Jfile))
            self.obfuscation.append(ObfLayer(int(inp), I, J))

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

    def obfuscate(self):
        self.obfuscation = [self._obfuscate_layer(layer) for layer in self.bp]

    def save(self, directory):
        assert self.obfuscation is not None
        if not os.path.exists(directory):
            os.mkdir(directory)
        self.ge.x0.save('%s/x0' % directory)
        self.ge.pzt.save('%s/pzt' % directory)
        for idx, layer in enumerate(self.obfuscation):
            Integer(layer.inp).save('%s/%d.input' % (directory, idx))
            layer.I.save('%s/%d.I' % (directory, idx))
            layer.J.save('%s/%d.J' % (directory, idx))

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
