#!/usr/bin/env sage -python

from __future__ import print_function
import unittest
import doctest

from sage.all import *

alpha = 16
beta = 8
eta = 138
kappa = 1
rho = 16
n = 138

ms = [long(1)] * n
    
def crt(a, b, m, n):
    g, alpha, beta = XGCD(m, n)
    q, r = Integer(b - a).quo_rem(g)
    if r != 0:
        raise ValueError("No solution to crt problem since gcd(%s,%s) does not divide %s-%s" % (
m, n, a, b))
    out = (a + q*alpha*m) % lcm(m, n)
    return out

def crt_list(v, moduli):
    x = v[0]
    m = moduli[0]
    for i in range(1,len(v)):
        x = crt(x,v[i],m,moduli[i])
        m = lcm(m,moduli[i])
        if i == 1:
            print("x[1] = %d" % x)
            print("m[1] = %d" % m)
    return x%m

import fastutils
x0, ps, gs, z, zinv, pzt \
  = fastutils.genparams(n, alpha, beta, eta, kappa)

min, max = 1 << rho - 1, (1 << rho) - 1
rs = [randint(min, max) for _ in xrange(n)]
rs = [long(r) for r in rs]
elems = [(r * g + m) * zinv % p for r, g, m, p in zip(rs, gs, ms, ps)]

r_py = crt_list(elems, ps)
print("*" * 90)
r_c = fastutils.encode(n, rho, ms, ps, gs, zinv, rs)
assert r_c == r_py
