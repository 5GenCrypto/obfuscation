#!/usr/bin/env sage -python

from sage.all import *
import sys


def logger(s, end='\n'):
    if end == '':
        print(s),
    else:
        print(s)
    sys.stdout.flush()


class GradedEncoding(object):
    def set_params(self, secparam, kappa):
        self.secparam = secparam
        self.kappa = kappa
        self.alpha = secparam
        self.beta = secparam
        self.rho = secparam
        self.mu = self.rho + self.alpha + self.secparam
        self.rho_f = self.kappa * (self.mu + self.rho + self.alpha + 2) + self.rho
        self.eta = self.rho_f + self.alpha + 2 * self.beta + self.secparam + 8
        self.nu = self.eta - self.beta - self.rho_f - self.secparam - 3
        self.n = self.eta * log(self.secparam, base=2)

    def print_params(self):
        print('Parameters:')
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

    def __init__(self, secparam=32, kappa=1):
        self.set_params(secparam, kappa)
        self.print_params()

        logger('Generating %d-bit primes p_i' % self.eta, end='')
        self.ps = []
        for _ in range(self.n):
            self.ps.append(random_prime((1 << self.eta) - 1))
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
        logger('Generating %d-bit primes g_i' % self.alpha, end='')
        self.gs = []
        for _ in range(self.n):
            self.gs.append(random_prime((1 << self.alpha) - 1))
            logger('.', end='')
        logger('')
        logger('Generating zero test element')
        zk = power_mod(self.z, self.kappa, self.x0)
        x0ps = [self.x0 / p for p in self.ps]
        gsis = [inverse_mod(g, p) for g, p in zip(self.gs, self.ps)]
        self.pzt = sum(randint(0, (1 << self.beta) - 1) * zk * gsis[i] * x0ps[i]
                       for i in range(self.n))

    def encode(self, msg):
        logger('Encoding value %d' % msg)
        ms = [0 for _ in range(self.n)]
        ms[0] = msg
        assert msg < self.gs[0], "Message must be smaller than g_0"
        logger('Generating random %d-bit integers r_i' % self.rho_f)
        rs = [random_prime((1 << self.rho_f) - 1) for _ in range(self.n)]
        logger('Generating elements for CRT')
        elems = [(r * g + m) * self.zinv % p
                 for r, g, m, p in zip(rs, self.gs, ms, self.ps)]
        logger('Finding c')
        return CRT(elems, self.ps)

    def is_zero(self, c):
        omega = (self.pzt * c) % self.x0
        return abs(omega) < (self.x0 >> self.nu)

    def add(self, cs):
        return reduce(operator.add, cs) % self.x0

    def mult(self, cs):
        return reduce(operator.mul, cs) % self.x0


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("invalid number of arguments")
        sys.exit(1)

    args = sys.argv[1].split()
    vals = []
    ops = []
    for arg in args:
        try:
            vals.append(int(arg))
        except ValueError:
            assert arg in ('+', '*'), "invalid argument"
            ops.append(arg)

    secparam = 16
    kappa = ops.count('*') + 1

    ge = GradedEncoding(secparam=secparam, kappa=kappa)
    encvals = [ge.encode(v) for v in vals]
    for op in ops:
        a, b = encvals.pop(0), encvals.pop(0)
        if op == '+':
            encvals.insert(0, ge.add([a, b]))
        else:
            encvals.insert(0, ge.mult([a, b]))
    assert len(encvals) == 1

    if ge.is_zero(encvals[0]):
        print('---- Output is ZERO')
    else:
        print('---- Output is something else')
