#!/usr/bin/env python2

from __future__ import print_function
from sage.all import *
import sys

SECPARAM = 4
KAPPA = 4

VECTOR_DIM = 20 * SECPARAM * KAPPA
ALPHA = SECPARAM
ETA = 10 * SECPARAM * KAPPA
RHO = SECPARAM


def log(s, **kwargs):
    print(s, **kwargs)
    sys.stdout.flush()


class GradedEncoding(object):
    def __init__(self):
        log('Generating %d %d-bit primes p_i' % (VECTOR_DIM, ETA), end='')
        self.ps = []
        for _ in range(VECTOR_DIM):
            self.ps.append(random_prime(2**ETA - 1))
            log('.', end='')
        log('')
        log('Computing x0', end='')
        self.x0 = reduce(operator.mul, self.ps)
        log('')
        log('Generating z', end='')
        while True:
            self.z = randint(0, self.x0)
            try:
                self.zinv = inverse_mod(self.z, self.x0)
                break
            except:
                pass
        log('')
        log('Generating %d %d-bit primes g_i' % (VECTOR_DIM, ALPHA), end='')
        self.gs = []
        for _ in range(VECTOR_DIM):
            self.gs.append(random_prime(2**ALPHA - 1))
            log('.', end='')
        log('')

    def encode(self, m):
        ms = [0 for _ in range(VECTOR_DIM)]
        ms[0] = m
        log('Generating %d random %d-bit integers r_i' % (VECTOR_DIM, RHO),
            end='')
        rs = []
        for _ in range(VECTOR_DIM):
            rs.append(random_prime(2**RHO - 1))
            log('.', end='')
        log('')
        log('Generating elements for CRT', end='')
        elems = []
        for i in range(VECTOR_DIM):
            elems.append(((rs[i] * self.gs[i] + ms[i]) * self.zinv) % self.ps[i])
        log('')
        log('Finding c', end='')
        c = CRT(elems, self.ps)
        log('')
        return c


if __name__ == '__main__':
    print('Parameters:')
    print('  Security Parameter: %d' % SECPARAM)
    print('  Kappa: %d' % KAPPA)
    print('  Vector Dimension: %d' % VECTOR_DIM)
    print('  Alpha: %d' % ALPHA)
    print('  Eta: %d' % ETA)
    print('  Rho: %d' % RHO)

    ge = GradedEncoding()
    c = ge.encode(0)
    print('encoded value = %d' % c)
