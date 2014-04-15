from __future__ import print_function

from branchingprogram import BranchingProgram, _group
import _obfuscator as _obf
import groups, utils

from sage.all import copy, VectorSpace, Zmod, ZZ

import collections, os, sys, time
import numpy as np

class ObfuscationException(Exception):
    pass

class AbstractObfuscator(object):
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
        if self._fast:
            self.n = self.eta
        else:
            self.n = int(self.eta * np.log2(self.secparam))
        self._print_params()

    def __init__(self, verbose=False, fast=False):
        self.obfuscation = None
        self._verbose = verbose
        self._fast = fast
        self.logger = utils.make_logger(self._verbose)
        if self._fast:
            self.logger('* Using small parameters for speed')

    def _gen_mlm_params(self, size, prime, nzs, directory):
        self.logger('Generating MLM parameters...')
        start = time.time()
        if not os.path.exists(directory):
            os.mkdir(directory)
        _obf.setup(size, self.n, self.alpha, self.beta, self.eta, self.nu,
                   self.rho, nzs, prime, directory)
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))

    def _obfuscate(self, bp):
        def _obfuscate_level(level, idx):
            self.logger('Obfuscating\n%s with set %s' % (
                level.zero, level.zeroset))
            self.logger('Obfuscating\n%s with set %s' % (
                level.one, level.oneset))
            start = time.time()
            # XXX: the C matrix multiplication code assumes the matrices are
            # stored column-wise, when in python they are stored row-wise, thus
            # the need for the transpose call.  This should be fixed.
            zero = [long(i) for i in level.zero.transpose().list()]
            one = [long(i) for i in level.one.transpose().list()]
            _obf.encode_level(idx, level.inp, zero, one, level.zeroset,
                              level.oneset)
            end = time.time()
            self.logger('Obfuscating level took: %f seconds' % (end - start))
        start = time.time()
        for idx, level in enumerate(bp):
            _obfuscate_level(level, idx)
        end = time.time()
        self.logger('Obfuscation took: %f seconds' % (end - start))

    def obfuscate(self, bp, secparam, directory, islayered):
        if bp.randomized:
            raise ObfuscationException('BPs must not be randomized')
        # remove old files in obfuscation directory
        if os.path.isdir(directory):
            files = os.listdir(directory)
            for file in os.listdir(directory):
                path = os.path.join(directory, file)
                os.unlink(path)
        # add two to kappa due to the bookend vectors
        kappa = len(bp) + 2
        self._set_params(secparam, kappa)
        prime = _obf.genprime(secparam)
        if islayered:
            bp.randomize(prime)
        else:
            R = Zmod(prime)
            alphas = [(R.random_element(), R.random_element())
                      for _ in xrange(len(bp))]
            bp.randomize(prime, alphas=alphas)
        # add two to nzs to take bookend vectors into account
        nzs = bp.set_straddling_sets() + 2
        self.logger('Number of Zs: %d' % nzs)
        if islayered:
            size = bp.size
        else:
            size = len(_group)
        self._gen_mlm_params(size, prime, nzs, directory)
        if not islayered:
            self._construct_multiplicative_constants(bp, alphas)
        self._construct_bookend_vectors(bp, prime, nzs)
        self._obfuscate(bp)

    def evaluate(self, directory, inp):
        start = time.time()
        files = os.listdir(directory)
        inputs = sorted(filter(lambda s: 'input' in s, files))
        islayered = True if type(self) == LayeredObfuscator else False
        result = _obf.evaluate(directory, inp, len(inputs), islayered)
        end = time.time()
        self.logger('Evaluation took: %f seconds' % (end - start))
        return result

class Obfuscator(AbstractObfuscator):
    def __init__(self, **kwargs):
        super(Obfuscator, self).__init__(**kwargs)

    def _construct_multiplicative_constants(self, bp, alphas):
        start = time.time()
        for idx, (layer, (a0, a1)) in enumerate(zip(bp, alphas)):
            _obf.encode_scalar(long(a0), layer.zeroset, "%d.a0_enc" % idx)
            _obf.encode_scalar(long(a1), layer.oneset, "%d.a1_enc" % idx)
        end = time.time()
        self.logger('Constructing multiplicative constants took: %f seconds' % (end - start))

    def _construct_bookend_vectors(self, bp, prime, nzs):
        start = time.time()
        sidx, tidx = nzs - 2, nzs - 1
        VSZp = VectorSpace(ZZ.residue_field(ZZ.ideal(prime)), _group.length)
        s = VSZp.random_element() * bp.m0i
        t = bp.m0 * VSZp.random_element()
        p = s * t
        _obf.encode_scalar(long(p), [sidx, tidx], "p_enc")
        _obf.encode_vector([long(i) for i in s], [sidx], "s_enc")
        _obf.encode_vector([long(i) for i in t], [tidx], "t_enc")
        end = time.time()
        self.logger('Constructing bookend vectors took: %f seconds' % (end - start))

    def obfuscate(self, bp, secparam, directory):
        super(Obfuscator, self).obfuscate(bp, secparam, directory, False)

class LayeredObfuscator(AbstractObfuscator):
    def __init__(self, **kwargs):
        super(LayeredObfuscator, self).__init__(**kwargs)

    def _construct_bookend_vectors(self, bp, prime, nzs):
        start = time.time()
        sidx, tidx = nzs - 2, nzs - 1
        s = bp.e_1 * bp.m0i
        t = bp.m0 * bp.e_w
        _obf.encode_vector([long(i) for i in s], [sidx], "s_enc")
        _obf.encode_vector([long(i) for i in t], [tidx], "t_enc")
        end = time.time()
        self.logger('Constructing bookend vectors took: %f seconds' % (end - start))

    def obfuscate(self, bp, secparam, directory):
        super(LayeredObfuscator, self).obfuscate(bp, secparam, directory, True)
