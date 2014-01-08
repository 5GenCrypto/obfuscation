#!/usr/bin/env python2
#
# Converts a circuit with AND and NOT gates to branching program.
#

from __future__ import print_function

import itertools
import numpy as np
import sys
from sage.all import *

MS = sage.matrix.matrix_space.MatrixSpace(GF(3), 3, 3)
MSZp = sage.matrix.matrix_space.MatrixSpace(
    ZZ.residue_field(ZZ.ideal(3388445611)), 3, 3)

I = MS.identity_matrix()
A = MS.matrix([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
Ai = A.inverse()
B = MS.matrix([[0, 0, 1], [0, -1, 0], [1, 0, 0]])
Bi = B.inverse()
C = MS.matrix([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
Ac = MS.matrix([[-1, 0, 0], [0, 0, -1], [0, -1, 0]])
Aci = Ac.inverse()
Bc = MS.matrix([[0, 1, 0], [1, 0, 1], [-1, 0, 1]])
Bci = Bc.inverse()
Cc = MS.matrix([[-1, 0, -1], [0, -1, 0], [0, 0, 1]])
Cci = Cc.inverse()


def ints(*args):
    return (int(arg) for arg in args)


def flatten(l):
    return list(itertools.chain(*l))


class Layer(object):
    def __init__(self, inp, I, J):
        self.inp = inp
        self.I = I
        self.J = J
    def __repr__(self):
        return "%d\nI:%s\nJ:%s" % (self.inp, self.I, self.J)


def prepend(layers, M):
    return [Layer(layers[0].inp, M * layers[0].I, M * layers[0].J)] \
        + layers[1:]
def append(layers, M):
    return layers[:-1] \
        + [Layer(layers[-1].inp, layers[-1].I * M, layers[-1].J * M)]


CONJUGATES = {
    'A': (Ac, Aci),
    'B': (Bc, Bci),
    'Ai': (Ac, Aci),
    'Bi': (Bc, Bci)
}


def conjugate(layers, target):
    if len(layers) == 1:
        l = layers[0]
        return [Layer(l.inp,
                      CONJUGATES[target][1] * l.I * CONJUGATES[target][0],
                      CONJUGATES[target][1] * l.J * CONJUGATES[target][0])]
    else:
        layers = prepend(layers, CONJUGATES[target][1])
        return append(layers, CONJUGATES[target][0])


def invert(layers):
    if len(layers) == 1:
        return [Layer(layers[0].inp,
                      np.linalg.inv(layers[0].I),
                      np.linalg.inv(layers[0].J))]
    else:
        return [invert(l) for l in reversed(layers)]


def notgate(layers):
    if len(layers) == 1:
        layer = layers[0]
        return [Layer(layer.inp,
                      Cci * layer.I * Cc * C,
                      Cci * layer.J * Cc * C)]
    else:
        layers = prepend(layers, Cci)
        return append(layers, Cc * C)


def circuit_to_bp(fname):
    ms = []
    with open(fname) as f:
        for line in f:
            if line.startswith('#'):
                continue
            num, rest = line.split(None, 1)
            num = int(num)
            if rest.startswith('input'):
                ms.append([Layer(num, I, C)])
            elif rest.startswith('gate') or rest.startswith('output'):
                # XXX: watchout!  I'm not sure what'll happen if we have
                # multiple outputs in a circuit
                if rest.startswith('gate'):
                    _, _, arity, _, rest = rest.split(None, 4)
                else:
                    _, _, _, arity, _, rest = rest.split(None, 5)
                arity = int(arity)
                if arity == 1:
                    _, a, b, _, _, _, in1, _ = rest.split(None, 7)
                    a, b, in1 = ints(a, b, in1)
                    if a == 1 and b == 0:
                        # NOT gate
                        a = notgate(ms[in1])
                        ms.append(a)
                    else:
                        raise("error: only support NOT so far:", line.strip())
                elif arity == 2:
                    _, a, b, c, d, _, _, _, in1, in2, _ \
                        = rest.split(None, 10)
                    a, b, c, d, in1, in2 = ints(a, b, c, d, in1, in2)
                    if a == b == c == 0 and d == 1:
                        # AND gate
                        a = conjugate(ms[in1], 'A')
                        b = conjugate(ms[in2], 'B')
                        c = conjugate(ms[in1], 'Ai')
                        d = conjugate(ms[in2], 'Bi')
                        ms.append(flatten([a, b, c, d]))
                    elif a == d == 0 and b == c == 1:
                        # XOR gate
                        ms.append(flatten([ms[in1], ms[in2]]))
                    else:
                        raise("error: unsupported gate:", line.strip())
                else:
                    raise("error: unsupported arity %d" % arity)
            else:
                raise("error: unknown type")
    return ms[-1]


def n_inputs(fname):
    nins = 0
    with open(fname) as f:
        for line in f:
            if line.startswith('#'):
                continue
            _, rest = line.split(None, 1)
            if rest.startswith('input'):
                nins = nins + 1
            else:
                break
    return nins


def obliviate(bp, nins, depth):
    newbp = []
    for m in bp:
        for i in range(nins):
            if m.inp == i:
                newbp.append(m)
            else:
                newbp.append(Layer(i, I, I))
    ms_needed = (4 ** depth) * nins
    for _ in range((ms_needed - len(newbp)) // nins):
        for i in range(nins):
            newbp.append(Layer(i, I, I))
    assert(len(newbp) == ms_needed)
    return newbp


def randomize(bp):
    def random_matrix():
        while True:
            m = MSZp.random_element()
            if not m.is_singular():
                return m, m.inverse()
    bp[0] = Layer(bp[0].inp,
                  MSZp.matrix(bp[0].I),
                  MSZp.matrix(bp[0].J))
    for i in range(1, len(bp)):
        mi, mii = random_matrix()
        bp[i-1] = Layer(bp[i-1].inp,
                        MSZp.matrix(bp[i-1].I) * mii,
                        MSZp.matrix(bp[i-1].J) * mii)
        bp[i] = Layer(bp[i].inp,
                      mi * MSZp.matrix(bp[i].I),
                      mi * MSZp.matrix(bp[i].J))
    bp[-1] = Layer(bp[-1].inp,
                   MSZp.matrix(bp[-1].I),
                   MSZp.matrix(bp[-1].J))
    return bp


def eval_bp(bp, inp, group):
    comp = group.identity_matrix()
    for m in bp:
        comp = comp * (m.I if inp[m.inp] == '0' else m.J)
    comp = MS.matrix(comp)
    if comp == I:
        return 0
    elif comp == C:
        return 1
    else:
        print("error: invalid return matrix:\n%s" % comp)
        return None

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("error: invalid arguments")
        exit(-1)

    fname = sys.argv[1]
    inp = sys.argv[2]
    bp = circuit_to_bp(fname)
    print("BRANCHING PROGRAM =", bp)
    print("BP length =", len(bp))
    # nins = n_inputs(fname)
    # bp = obliviate(bp, nins, 1)
    bp = randomize(bp)
    print("RANDOMIZED BP =", bp)
    print("RBP length =", len(bp))
    out = eval_bp(bp, inp, MSZp)
    print("OUTPUT = %d" % out)
