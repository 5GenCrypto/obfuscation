#!/usr/bin/env sage -python

from __future__ import print_function

from sage.all import *

import math, sys, time
import functools
import utils

def genprimes(num, bitlength, ncpus):
    @parallel(ncpus=ncpus)
    def _random_prime(r):
        with seed(r):
            return random_prime(bitlength)
    rs = [randint(0, 1 << 64) for _ in xrange(num)]
    ps = list(_random_prime(rs))
    return [p for _, p in ps]

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
        # FIXME: we're currently using a too small secparam, and thus we need
        # alpha = 2 * secparam otherwise encoding fails due to message being
        # smaller than the g_i values
        self.alpha = secparam << 1
        self.beta = secparam
        self.rho = 2 * secparam
        self.mu = self.rho + self.alpha + self.secparam
        self.rho_f = self.kappa * (self.mu + self.rho + self.alpha + 2) + self.rho
        self.eta = self.rho_f + self.alpha + 2 * self.beta + self.secparam + 8
        self.nu = self.eta - self.beta - self.rho_f - self.secparam - 3
        # XXX: things appear to break when we use the real value for n.
        # However, smaller values (e.g., n = eta) appear to work...
        self.n = self.eta * 2
        # self.n = int(self.eta * math.log(self.secparam, 2))

    def _print_params(self):
        print('Graded Encoding Parameters:')
        print('  Lambda: %d' % self.secparam)
        print('  Kappa: %d' % self.kappa)
        print('  Alpha: %d' % self.alpha)
        print('  Beta: %d' % self.beta)
        print('  Eta: %d' % self.eta)
        print('  Mu: %d' % self.mu)
        print('  Nu: %d' % self.nu)
        print('  Rho: %d' % self.rho)
        print('  Rhof: %d' % self.rho_f)
        print('  N: %d' % self.n)

    def __init__(self, verbose=False, parallel=False, ncpus=1):
        self._verbose = verbose
        self._parallel = parallel
        self._ncpus = ncpus
        self.logger = functools.partial(utils.logger, verbose=self._verbose)

    def load_system_params(self, secparam, kappa, x0, pzt):
        self._set_params(secparam, kappa)
        if self._verbose:
            self._print_params()
        self.x0 = x0
        self.pzt = pzt

    def gen_system_params(self, secparam, kappa):
        self._set_params(secparam, kappa)
        if self._verbose:
            self._print_params()

        self.logger('Generating %d-bit primes p_i ' % self.eta, end='')
        start = time.time()
        primesize = (1 << self.eta) - 1
        if self._parallel:
            self.logger('(parallel)')
            self.ps = genprimes(self.n, primesize, self._ncpus)
        else:
            self.logger('(sequential)')
            self.ps = [random_prime(primesize) for _ in xrange(self.n)]
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))

        self.logger('Computing x0')
        start = time.time()
        self.x0 = reduce(operator.mul, self.ps)
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))

        self.logger('Generating z')
        start = time.time()
        self.z, self.zinv = genz(self.x0)
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))

        self.logger('Generating %d-bit primes g_i (sequentially)' % self.alpha)
        start = time.time()
        primesize = (1 << self.alpha) - 1
        self.gs = [random_prime(primesize) for _ in xrange(self.n)]
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))

        self.logger('Generating zero test element')
        start = time.time()
        zk = power_mod(self.z, self.kappa, self.x0)
        end = time.time()
        self.logger('  Computing power_mod: %f seconds' % (end - start))
        start = time.time()
        x0ps = [self.x0 / p for p in self.ps]
        end = time.time()
        self.logger('  Computing x0 / p_i: %f seconds' % (end - start))
        start = time.time()
        gsis = [inverse_mod(g, p) for g, p in zip(self.gs, self.ps)]
        end = time.time()
        self.logger('  Computing inverse_mod: %f seconds' % (end - start))
        start = time.time()
        self.pzt = sum(randint(0, (1 << self.beta) - 1) * zk * gsis[i] * x0ps[i]
                       for i in xrange(self.n))
        end = time.time()
        self.logger('  Computing pzt: %f seconds' % (end - start))

    def encode(self, val):
        self.logger('Encoding value %d' % val)
        ms = [0] * self.n
        ms[0] = val
        assert val < self.gs[0], "Message must be smaller than g_0"
        min, max = 1 << self.rho - 1, (1 << self.rho) - 1
        start = time.time()
        rs = [randint(min, max) for _ in xrange(self.n)]
        end = time.time()
        self.logger('  Generating r_i values: %f seconds' % (end - start))
        start = time.time()
        elems = [(r * g + m) * self.zinv % p
                 for r, g, m, p in zip(rs, self.gs, ms, self.ps)]
        end = time.time()
        self.logger('  Generating elements: %f seconds' % (end - start))
        start = time.time()
        r = CRT(elems, self.ps)
        end = time.time()
        self.logger('  CRT: %f seconds' % (end - start))
        return r

    def encode_list(self, vals):
        @parallel(ncpus=self._ncpus)
        def _encode(val, r):
            with seed(r):
                return self.encode(val)
        self.logger('Encoding values %s' % vals)
        rs = [randint(0, 1 << 64) for _ in xrange(len(vals))]
        es = list(_encode(zip(vals, rs)))
        return [e for _, e in es]

    def is_zero(self, c):
        omega = (self.pzt * c) % self.x0
        return abs(omega) < (self.x0 >> self.nu)

    def add(self, cs):
        return reduce(operator.add, cs) % self.x0

    # def sub(self, cs):
    #     return reduce(operator.sub, cs) % self.x0

    def mult(self, cs):
        return reduce(operator.mul, cs) % self.x0
