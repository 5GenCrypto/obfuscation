#!/usr/bin/env sage -python

from __future__ import print_function
from sage.all import *
import time

N = 1064
B = 266
LENGTH = (1 << B) - 1

print('Generating primes (sequentially)')
start = time.time()
rs = [randint(0, 1 << 64) for _ in xrange(N)]
ps = []
for r in rs:
    with seed(r):
        ps.append(random_prime(LENGTH))
end = time.time()
print('Took: %f seconds' % (end - start))

print('Number of CPUs = %d' % sage.parallel.ncpus.ncpus())

@parallel(ncpus=sage.parallel.ncpus.ncpus())
def _random_prime(r):
    with seed(r):
        return random_prime(LENGTH)
print('Generating primes (in parallel)')
start = time.time()
rs = [randint(0, 1 << 64) for _ in xrange(N)]
ps = list(_random_prime(rs))
end = time.time()
print('Took: %f seconds' % (end - start))
ps = [p for _, p in ps]
