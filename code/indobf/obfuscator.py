#!/usr/bin/env sage -python

from __future__ import print_function

from branchingprogram import BranchingProgram
import utils, fastutils

from sage.all import (flatten, Integer, load, MatrixSpace, randint,
                      random_prime, vector, VectorSpace, Zmod, ZZ)

import collections, os, sys, time
import numpy

def ms2list(m):
    '''Convert an element in MS to a flat integer list'''
    m = [[long(e) for e in row] for row in m]
    return [long(e) for e in flatten(m)]

def to_long(l):
    return [long(e) for e in l]

ObfLayer = collections.namedtuple('ObfLayer', ['inp', 'zero', 'one'])

def load_layer(directory, inp, zero, one):
    inp = load('%s/%s' % (directory, inp))
    zero = load('%s/%s' % (directory, zero))
    one = load('%s/%s' % (directory, one))
    return ObfLayer(int(inp), zero, one)

def save_layer(layer, directory, idx):
    Integer(layer.inp).save('%s/%d.input' % (directory, idx))
    layer.zero.save('%s/%d.zero' % (directory, idx))
    layer.one.save('%s/%d.one' % (directory, idx))

class ObfuscationException(Exception):
    pass

class Obfuscator(object):

    def _print_params(self):
        self.logger('Graded encoding parameters:')
        self.logger('  Security Parameter: %d' % self.secparam)
        self.logger('  Kappa: %d' % self.kappa)
        self.logger('  Alpha: %d' % self.alpha)
        self.logger('  Beta: %d' % self.beta)
        self.logger('  Eta: %d' % self.eta)
        self.logger('  Nu: %d' % self.nu)
        self.logger('  Rho: %d' % self.rho)
        self.logger('  Rho_f: %d' % self.rho_f)
        self.logger('  N: %d' % self.n)

    def _set_params(self, secparam, kappa):
        self.secparam = secparam
        self.kappa = kappa
        self.alpha = self.secparam
        self.beta = self.secparam
        self.rho = self.secparam
        self.rho_f = self.kappa * (self.rho + self.alpha + 2)
        self.eta = self.rho_f + self.alpha + 2 * self.beta + self.secparam + 8
        self.nu = self.eta - self.beta - self.rho_f - self.secparam - 3
        assert self.nu >= self.alpha + self.beta + 5
        # XXX: use smaller n value for now to speed things up
        if self._fast:
            self.n = self.eta
        else:
            self.n = int(self.eta * numpy.log2(self.secparam))
        self._print_params()

    def __init__(self, verbose=False, disable_mbundling=False,
                 disable_bookends=False, fast=False):
        self.obfuscation = None
        self._verbose = verbose
        self._disable_mbundling = disable_mbundling
        self._disable_bookends = disable_bookends
        self._fast = fast

        self.logger = utils.make_logger(self._verbose)
        if self._disable_mbundling:
            self.logger('* Multiplicative bundling disabled')
        if self._disable_bookends:
            self.logger('* Bookends disabled')
        if self._fast:
            self.logger('* Using small parameters for speed')

    def save(self, directory):
        assert self.obfuscation is not None

        if not os.path.exists(directory):
            os.mkdir(directory)
        Integer(self.nu).save('%s/nu' % directory)
        Integer(self.x0).save('%s/x0' % directory)
        Integer(self.pzt).save('%s/pzt' % directory)
        if not self._disable_bookends:
            Integer(self.p_enc).save('%s/p_enc' % directory)
            vector(self.s_enc).save('%s/s_enc' % directory)
            vector(self.t_enc).save('%s/t_enc' % directory)
        if not self._disable_mbundling:
            vector(self.a0s_enc).save('%s/a0s_enc' % directory)
            vector(self.a1s_enc).save('%s/a1s_enc' % directory)
        for idx, layer in enumerate(self.obfuscation):
            save_layer(layer, directory, idx)

    def load(self, directory):
        assert self.obfuscation is None

        self.nu = int(load('%s/nu.sobj' % directory))
        self.x0 = long(load('%s/x0.sobj' % directory))
        pzt = long(load('%s/pzt.sobj' % directory))
        if not self._disable_bookends:
            self.p_enc = long(load('%s/p_enc.sobj' % directory))
            self.s_enc = load('%s/s_enc.sobj' % directory)
            self.s_enc = [long(e) for e in self.s_enc]
            self.t_enc = load('%s/t_enc.sobj' % directory)
            self.t_enc = [long(e) for e in self.t_enc]
        if not self._disable_mbundling:
            self.a0s_enc = load('%s/a0s_enc.sobj' % directory)
            self.a0s_enc = [long(e) for e in self.a0s_enc]
            self.a1s_enc = load('%s/a1s_enc.sobj' % directory)
            self.a1s_enc = [long(e) for e in self.a1s_enc]
        files = os.listdir(directory)
        inputs = sorted(filter(lambda s: 'input' in s, files))
        zeros = sorted(filter(lambda s: 'zero' in s, files))
        ones = sorted(filter(lambda s: 'one' in s, files))
        self.obfuscation = [load_layer(directory, inp, zero, one) for inp, zero,
                            one in zip(inputs, zeros, ones)]
        fastutils.loadparams(self.x0, pzt)

    def _gen_mlm_params(self, prime, nzs):
        self.logger('Generating MLM parameters...')
        start = time.time()
        self.x0, self.pzt = fastutils.genparams(self.n, self.alpha, self.beta,
                                                self.eta, self.kappa, self.rho,
                                                nzs, [prime])
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))

    def _construct_multiplicative_constants(self, bp, alphas):
        if not self._disable_mbundling:
            self.logger('Constructing multiplicative constants...')
            start = time.time()
            self.a0s_enc = []
            self.a1s_enc = []
            for layer, (a0, a1) in zip(bp, alphas):
                a0enc = fastutils.encode_scalar(long(a0), 0, layer.zeroset)
                a1enc = fastutils.encode_scalar(long(a1), 0, layer.oneset)
                self.a0s_enc.append(a0enc)
                self.a1s_enc.append(a1enc)
            end = time.time()
            self.logger('Took: %f seconds' % (end - start))

    def _construct_bookend_vectors(self, bp, prime, nzs, veclen):
        if not self._disable_bookends:
            self.logger('Constructing bookend vectors...')
            start = time.time()
            sidx, tidx = nzs - 2, nzs - 1
            VSZp = VectorSpace(ZZ.residue_field(ZZ.ideal(prime)), veclen)
            s = VSZp.random_element() * bp.m0i
            t = bp.m0 * VSZp.random_element()
            p = s * t
            self.p_enc = fastutils.encode_scalar(long(p), 0, [sidx, tidx])
            if self._disable_mbundling:
                for i in xrange(nzs - 2):
                    self.p_enc *= fastutils.encode_scalar(1L, 0, [i])
            self.s_enc = fastutils.encode_vector([long(i) for i in s], 0, [sidx])
            self.t_enc = fastutils.encode_vector([long(i) for i in t], 0, [tidx])
            end = time.time()
            self.logger('Took: %f seconds' % (end - start))

    def _obfuscate(self, bp, veclen):
        length = veclen * veclen * 2
        half = length >> 1
        def _obfuscate_layer(layer):
            self.logger('Obfuscating\n%s with set %s' % (
                layer.zero, layer.zeroset))
            self.logger('Obfuscating\n%s with set %s' % (
                    layer.one, layer.oneset))
            start = time.time()
            m = ms2list(layer.zero)
            m.extend(ms2list(layer.one))
            es = fastutils.encode_layer(m, 0, layer.zeroset, layer.oneset)
            zero, one = self.MS(es[:half]), self.MS(es[half:])
            end = time.time()
            self.logger('Obfuscating layer took: %f seconds' % (end - start))
            return ObfLayer(layer.inp, zero, one)
        self.logger('Obfuscating...')
        start = time.time()
        self.obfuscation = [_obfuscate_layer(layer) for layer in bp]
        end = time.time()
        self.logger('Obfuscation took: %f seconds' % (end - start))

    def obfuscate(self, bp, secparam):
        group = bp.bpgroup
        if bp.randomized:
            raise ObfuscationException('BPs must not be randomized')
        
        self.MS = MatrixSpace(ZZ, group.length)

        kappa = len(bp)
        if not self._disable_bookends:
            # add two to kappa due to the bookend vectors
            kappa += 2

        self._set_params(secparam, kappa)

        prime = long(random_prime((1 << secparam) - 1,
                                  lbound=(1 << secparam - 1)))

        if self._disable_mbundling:
            alphas = None
        else:
            R = Zmod(prime)
            alphas = [(R.random_element(), R.random_element())
                      for _ in xrange(len(bp))]
        bp.randomize(prime, alphas=alphas)

        nzs = bp.set_straddling_sets()
        if not self._disable_bookends:
            # take bookend vectors into account
            nzs = nzs + 2
        self.logger('Number of Zs: %d' % nzs)

        self._gen_mlm_params(prime, nzs)

        self._construct_multiplicative_constants(bp, alphas)

        self._construct_bookend_vectors(bp, prime, nzs, group.length)

        self._obfuscate(bp, group.length)

    def _is_zero(self, c):
        return fastutils.is_zero(long(c), self.nu)
        
    def evaluate(self, inp):
        assert self.obfuscation is not None

        start = time.time()

        first = self.obfuscation[0]
        p1 = first.zero if inp[first.inp] == '0' else first.one
        for m in self.obfuscation[1:]:
            p1 *= m.zero if inp[m.inp] == '0' else m.one
        if self._disable_bookends:
            result = 0 if self._is_zero(p1[0][1]) and self._is_zero(p1[1][0]) else 1
        else:
            # need to use numpy arrays here, as sage constructs cause weird issues
            p1 = long((numpy.dot(numpy.dot(self.s_enc, numpy.array(p1)),
                                 self.t_enc)) % self.x0)

            p2 = self.p_enc
            if not self._disable_mbundling:
                for i, m in enumerate(self.obfuscation):
                    p2 = p2 * (self.a0s_enc[i] if inp[m.inp] == '0' else self.a1s_enc[i])
        
            result = 0 if self._is_zero(p1 - p2) else 1

        end = time.time()

        self.logger('Evaluation took: %f seconds' % (end - start))

        return result
