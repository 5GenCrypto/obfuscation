#!/usr/bin/env sage -python

from __future__ import print_function

from sage.all import *

import math, sys, time
import utils
import fastutils

def genz(x0):
    while True:
        z = randint(0, x0)
        try:
            zinv = inverse_mod(z, x0)
            return z, zinv
        except:
            pass

class GradedEncoding(object):

    def _set_params(self, secparam, kappa):
        self.secparam = secparam
        self.kappa = kappa
        if secparam <= 16:
            # for small security parameters we need alpha = 2 * secparam,
            # otherwise encoding fails due to message being smaller than the g_i
            # values
            self.alpha = secparam << 1
        else:
            self.alpha = secparam
        self.beta = secparam
        self.rho = 2 * secparam
        mu = self.rho + self.alpha + self.secparam
        rho_f = self.kappa * (mu + self.rho + self.alpha + 2) + self.rho
        # rho_f = self.kappa * (self.rho + self.alpha + 2) + self.rho
        self.eta = rho_f + self.alpha + 2 * self.beta + self.secparam + 8
        self.nu = self.eta - self.beta - rho_f - self.secparam - 3
        # XXX: use smaller n value for now to speed things up
        self.n = self.eta
        # self.n = int(self.eta * math.log(self.secparam, 2))

    def _print_params(self):
        print('Graded encoding parameters:')
        print('  Lambda: %d' % self.secparam)
        print('  Kappa: %d' % self.kappa)
        print('  Alpha: %d' % self.alpha)
        print('  Beta: %d' % self.beta)
        print('  Eta: %d' % self.eta)
        print('  Nu: %d' % self.nu)
        print('  Rho: %d' % self.rho)
        print('  N: %d' % self.n)

    def genprimes(self, num, bitlength):
        @parallel(ncpus=self._ncpus)
        def _random_prime(r):
            with seed(r):
                return random_prime(bitlength, proof=False)
        if self._parallel:
            rs = [randint(0, 1 << 64) for _ in xrange(num)]
            ps = list(_random_prime(rs))
            return [p for _, p in ps]
        else:
            return [random_prime(bitlength) for _ in xrange(num)]

    def __init__(self, verbose=False, parallel=False, ncpus=1, use_c=False):
        self._verbose = verbose
        self._parallel = parallel
        self._ncpus = ncpus
        self._use_c = use_c
        self.logger = utils.make_logger(self._verbose)

    def load_system_params(self, secparam, kappa, x0, pzt):
        self._set_params(secparam, kappa)
        if self._verbose:
            self._print_params()
        self.x0 = x0
        self.pzt = pzt
        if self._use_c:
            fastutils.loadparams(long(x0), long(pzt))

    def gen_system_params(self, secparam, kappa):
        self._set_params(secparam, kappa)
        if self._verbose:
            self._print_params()

        self.logger('Generating parameters...')
        start = time.time()
        if self._use_c:
            self.x0, self.pzt = fastutils.genparams(self.n, self.alpha,
                                                    self.beta, self.eta,
                                                    self.kappa)
        else:
            primesize = (1 << self.eta) - 1
            self.ps = self.genprimes(self.n, primesize)
            self.x0 = reduce(operator.mul, self.ps)
            primesize = (1 << self.alpha) - 1
            self.gs = self.genprimes(self.n, primesize)
            self.z, self.zinv = genz(self.x0)
            zk = power_mod(self.z, self.kappa, self.x0)
            x0ps = [self.x0 / p for p in self.ps]
            gsis = [inverse_mod(g, p) for g, p in zip(self.gs, self.ps)]
            self.pzt = sum(randint(0, (1 << self.beta) - 1) * zk * gsis[i] * x0ps[i]
                           for i in xrange(self.n))
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))

    def encode(self, val):
        self.logger('Encoding value %d' % val)
        start = time.time()
        if self._use_c:
            r = fastutils.encode(long(val), self.rho)
        else:
            ms = [0] * self.n
            ms[0] = val
            assert ms[0] < self.gs[0], "Message must be smaller than g_0"
            min, max = 1 << self.rho - 1, (1 << self.rho) - 1
            rs = [randint(min, max) for _ in xrange(self.n)]
            elems = [(r * g + m) * self.zinv % p
                     for r, g, m, p in zip(rs, self.gs, ms, self.ps)]
            r = CRT_list(elems, self.ps)
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))
        return r

    def encode_list(self, vals):
        if self._use_c:
            vals = [long(val) for val in vals]
            return fastutils.encode_list(vals, self.rho)
        else:
            # _encode takes an index value as its first argument so we can recreate
            # the proper ordering of the encoded values after parallelization
            @parallel(ncpus=self._ncpus)
            def _encode(idx, val, r):
                with seed(r):
                    return idx, self.encode(val)
            self.logger('Encoding values %s' % vals)
            indices = list(xrange(len(vals)))
            # generate random values needed for encoding
            rs = [randint(0, 1 << 64) for _ in xrange(len(vals))]
            # compute encoded values in parallel
            es = list(_encode(zip(indices, vals, rs)))
            # extract encoded values from parallelization output
            es = [e for _, e in es]
            # sort by indices
            es.sort()
            # return encoded values, ignoring indices
            return [e for _, e in es]

    def is_zero(self, c):
        if self._use_c:
            return fastutils.is_zero(long(c), self.nu)
        else:
            omega = (self.pzt * c) % self.x0
            return abs(omega) < (self.x0 >> self.nu)

    def add(self, cs):
        return reduce(operator.add, cs) % self.x0

    def mult(self, cs):
        return reduce(operator.mul, cs) % self.x0
