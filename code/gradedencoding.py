#!/usr/bin/env python2

from __future__ import print_function
from sage.all import *
import sys

LAMBDA = 16
KAPPA = 1

ALPHA = LAMBDA
BETA = LAMBDA
RHO = LAMBDA
MU = RHO + ALPHA + LAMBDA
RHOf = KAPPA * (MU + RHO + ALPHA + 2) + RHO
ETA = RHOf + ALPHA + 2 * BETA + LAMBDA + 8
NU = ETA - BETA - RHOf - LAMBDA - 3
N = ETA * log(LAMBDA, base=2)


def print_params():
    print('Parameters:')
    print('  Lambda: %d' % LAMBDA)
    print('  Kappa: %d' % KAPPA)
    print('  Alpha: %d' % ALPHA)
    print('  Beta: %d' % BETA)
    print('  Eta: %d' % ETA)
    print('  Mu: %d' % MU)
    print('  Nu: %d' % NU)
    print('  Rho: %d' % RHO)
    print('  Rhof: %d' % RHOf)
    print('  N: %d' % N)


def logger(s, **kwargs):
    print(s, **kwargs)
    sys.stdout.flush()


class GradedEncoding(object):
    def __init__(self):
        logger('Generating %d %d-bit primes p_i' % (N, ETA), end='')
        self.ps = []
        for _ in range(N):
            self.ps.append(random_prime((1 << ETA) - 1))
            logger('.', end='')
        logger('')
        logger('Computing x0', end='')
        self.x0 = reduce(operator.mul, self.ps)
        logger('')
        logger('Generating z', end='')
        while True:
            self.z = randint(0, self.x0)
            try:
                self.zinv = inverse_mod(self.z, self.x0)
                break
            except:
                pass
        logger('')
        logger('Generating %d %d-bit primes g_i' % (N, ALPHA), end='')
        self.gs = []
        for _ in range(N):
            self.gs.append(random_prime((1 << ALPHA) - 1))
            logger('.', end='')
        logger('')
        logger('Generating H', end='')
        nrows, ncols = floor(N / 2.0), ceil(N / 2.0)
        MS = MatrixSpace(ZZ, nrows, ncols)
        Ifloor = MS(identity_matrix(nrows))
        Iceil = MS(identity_matrix(ncols))
        Z = MS(zero_matrix(nrows, ncols))
        b = max((floor(BETA / ceil(log(1 + ncols, base=2))), 1))
        Hs = []
        for _ in range(b):
            A = MS(random_matrix(Zmod(3), nrows, ncols)) - \
                ones_matrix(nrows, ncols)
            HA = block_matrix(ZZ, 2, 2, [Ifloor, A, Iceil, Z])
            if randint(0, 1) == 1:
                HA = HA.transpose()
            Hs.append(HA)
        H = reduce(operator.mul, Hs)
        logger('')
        logger('Generating zero test vector', end='')
        self.pzts = []
        zk = power_mod(self.z, KAPPA, self.x0)
        x0ps = [self.x0 / p for p in self.ps]
        gsis = [inverse_mod(self.gs[i], self.ps[i])
                for i in range(len(self.gs))]
        for j in range(N):
            v = sum(H[i][j] * zk * gsis[i] * x0ps[i] for i in range(N))
            self.pzts.append(v % self.x0)
            logger('.', end='')
        logger('')

    def encode(self, m):
        ms = [0 for _ in range(N)]
        ms[0] = m
        logger('Generating %d random %d-bit integers r_i' % (N, RHO), end='')
        rs = []
        for _ in range(N):
            rs.append(random_prime((1 << RHO) - 1))
            logger('.', end='')
        logger('')
        logger('Generating elements for CRT', end='')
        elems = [(rs[i] * self.gs[i] + ms[i]) * self.zinv % self.ps[i]
                 for i in range(len(ms))]
        logger('')
        logger('Finding c', end='')
        c = CRT(elems, self.ps)
        logger('')
        return c

    def is_zero(self, c):
        omega = [c * e % self.x0 for e in self.pzts]
        return max(omega) < (self.x0 >> NU)


if __name__ == '__main__':
    print_params()

    if len(sys.argv) == 1:
        inp = 0
    else:
        inp = int(sys.argv[1])

    ge = GradedEncoding()
    print('---- Encoding value %d' % inp)
    c = ge.encode(inp)
    print('---- Encoded value = %d' % c)
    if ge.is_zero(c):
        print('---- Output is *0*')
    else:
        print('---- Output is *1*')
