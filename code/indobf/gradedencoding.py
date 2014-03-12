#!/usr/bin/env sage -python

from __future__ import print_function

from sage.all import *

import math, sys, time
import utils
import fastutils

class GradedEncoding(object):

    def _set_params(self, secparam, kappa):
        self.secparam = secparam
        self.kappa = kappa
        self.alpha = secparam
        self.beta = secparam
        self.rho = secparam
        # mu = self.rho + self.alpha + self.secparam
        # self.rho_f = self.kappa * (mu + self.rho + self.alpha + 2) + self.rho
        self.rho_f = self.kappa * (self.rho + self.alpha)
        self.eta = self.rho_f + self.alpha + 2 * self.beta + self.secparam + 8
        self.nu = self.eta - self.beta - self.rho_f - self.secparam - 3
        assert self.nu >= self.alpha + self.beta + 5
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
        print('  Rho_f: %d' % self.rho_f)
        print('  N: %d' % self.n)

    def __init__(self, verbose=False):
        self._verbose = verbose
        self.logger = utils.make_logger(self._verbose)

    def load_system_params(self, secparam, kappa, x0, pzt):
        self._set_params(secparam, kappa)
        if self._verbose:
            self._print_params()
        fastutils.loadparams(long(x0), long(pzt))

    def gen_system_params(self, secparam, kappa):
        self._set_params(secparam, kappa)
        if self._verbose:
            self._print_params()

        self.logger('Generating parameters...')
        start = time.time()
        self.x0, self.pzt = fastutils.genparams(self.n, self.alpha, self.beta,
                                                self.eta, self.kappa)
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))

    def encode(self, val):
        self.logger('Encoding value %d' % val)
        start = time.time()
        r = fastutils.encode(long(val), self.rho)
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))
        return r

    def encode_list(self, vals):
        vals = [long(val) for val in vals]
        return fastutils.encode_list(vals, self.rho)

    def is_zero(self, c):
        return fastutils.is_zero(long(c), self.nu)
