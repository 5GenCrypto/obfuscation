#!/usr/bin/env sage -python

from sage.all import *
import sys

LAMBDA = 16
KAPPA = 2

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


def logger(s, end='\n'):
    if end == '':
        print(s),
    else:
        print(s)
    sys.stdout.flush()


class GradedEncoding(object):
    def __init__(self):
        logger('Generating %d %d-bit primes p_i' % (N, ETA), end='')
        self.ps = []
        for _ in range(N):
            self.ps.append(random_prime((1 << ETA) - 1))
            logger('.', end='')
        logger('')
        logger('Computing x0')
        self.x0 = reduce(operator.mul, self.ps)
        logger('Generating z')
        while True:
            self.z = randint(0, self.x0)
            try:
                self.zinv = inverse_mod(self.z, self.x0)
                break
            except:
                pass
        logger('Generating %d %d-bit primes g_i' % (N, ALPHA), end='')
        self.gs = []
        for _ in range(N):
            self.gs.append(random_prime((1 << ALPHA) - 1))
            logger('.', end='')
        logger('')
        # logger('Generating H')
        # nrows, ncols = floor(N / 2.0), ceil(N / 2.0)
        # MS = MatrixSpace(ZZ, nrows, ncols)
        # Ifloor = MS(identity_matrix(nrows))
        # Iceil = MS(identity_matrix(ncols))
        # Z = MS(zero_matrix(nrows, ncols))
        # b = max((floor(BETA / ceil(log(1 + ncols, base=2))), 1))
        # Hs = []
        # for _ in range(b):
        #     A = MS(random_matrix(Zmod(3), nrows, ncols)) - \
        #         ones_matrix(nrows, ncols)
        #     HA = block_matrix(ZZ, 2, 2, [Ifloor, A, Iceil, Z])
        #     if randint(0, 1) == 1:
        #         HA = HA.transpose()
        #     Hs.append(HA)
        # H = matrix(reduce(operator.mul, Hs))
        # print('Checking singularity...')
        # assert not H.is_singular(), "H is singular!"
        # print('Checking norms...')
        # assert numpy.linalg.norm(H, ord=numpy.inf) <= (1 << BETA)
        # assert numpy.linalg.norm(H.inverse(), ord=numpy.inf) <= (1 << BETA)
        logger('Generating zero test element')
        zk = power_mod(self.z, KAPPA, self.x0)
        x0ps = [self.x0 / p for p in self.ps]
        gsis = [inverse_mod(g, p) for g, p in zip(self.gs, self.ps)]
        self.pzt = sum(randint(0, (1 << BETA) - 1) * zk * gsis[i] * x0ps[i]
                       for i in range(N))
        # logger('Generating zero test vector', end='')
        # self.pzts = []
        # for j in range(N):
        #     v = sum(H[i][j] * zk * gsis[i] * x0ps[i] for i in range(N))
        #     self.pzts.append(v % self.x0)
        #     logger('.', end='')
        # logger('')

    def encode(self, msg):
        ms = [0 for _ in range(N)]
        ms[0] = msg
        assert msg < self.gs[0], "Message must be smaller than g_0"
        logger('Generating %d random %d-bit integers r_i' % (N, RHOf))
        rs = [random_prime((1 << RHOf) - 1) for _ in range(N)]
        logger('Generating elements for CRT')
        elems = [(r * g + m) * self.zinv % p
                 for r, g, m, p in zip(rs, self.gs, ms, self.ps)]
        logger('Finding c')
        return CRT(elems, self.ps)

    def is_zero(self, c):
        omega = self.pzt * c % self.x0
        return abs(omega) < (self.x0 >> NU)
        # omega = [c * e % self.x0 for e in self.pzts]
        # return max(omega) < (self.x0 >> NU)


if __name__ == '__main__':
    print_params()

    if len(sys.argv) == 1:
        inp = [0]
    else:
        inp = [int(i) for i in sys.argv[1:]]

    ge = GradedEncoding()
    cs = []
    for i in inp:
        print('---- Encoding value %d' % i)
        c = ge.encode(i)
        cs.append(c)
    r = reduce(operator.mul, cs) % ge.x0
    if ge.is_zero(r):
        print('---- Output is ZERO')
    else:
        print('---- Output is something else')
