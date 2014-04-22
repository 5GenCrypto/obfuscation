#!/usr/bin/env sage -python

from __future__ import print_function
import functools, math, sys
import numpy as np

def logger(s, end='\n', verbose=False):
    if verbose:
        print(s, end=end, file=sys.stderr)
        sys.stderr.flush()

def make_logger(verbose):
    return functools.partial(logger, verbose=verbose)

#
# Below code taken from:
# https://stackoverflow.com/questions/4287721/easiest-way-to-perform-modular-matrix-inversion-with-python
#
# XXX: I haven't verified the correctness of this code!
#

# Finds the inverse of matrix m mod p
def mod_mat_inv(m, p):
  n = len(m)
  m = np.matrix(m)
  adj = np.zeros(shape=(n, n))
  for i in xrange(n):
    for j in xrange(n):
      adj[i][j] = ((-1)**(i+j) * int(round(np.linalg.det(minor(m, j, i))))) % p
  return (mod_inv(int(round(np.linalg.det(m))), p) * adj) % p

# Finds the inverse of a mod p, if it exists
def mod_inv(a, p):
  for i in range(1, p):
    if (i * a) % p == 1:
      return i
  raise ValueError(str(a) + " has no inverse mod " + str(p))

# Return matrix m with the ith row and jth column deleted
def minor(m, i, j):
  m = np.array(m)
  minor = np.zeros(shape=(len(m)-1, len(m)-1))
  p = 0
  for s in xrange(len(minor)):
    if p == i:
      p += 1
    q = 0
    for t in xrange(len(minor)):
      if q == j:
        q += 1
      minor[s][t] = m[p][q]
      q += 1
    p += 1
  return minor
