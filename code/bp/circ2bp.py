#!/usr/bin/env python3
#
# Converts a circuit with AND and NOT gates to branching program.
#

import itertools
import numpy as np
import sys

I = np.eye(3)
A = np.matrix('-1 0 0; 0 -1 0; 0 0 1')
Ai = np.linalg.inv(A)
B = np.matrix('0 0 1; 0 -1 0; 1 0 0')
Bi = np.linalg.inv(B)
C = np.matrix('-1 0 0; 0 1 0; 0 0 -1')
Ac = np.matrix('-1 0 0; 0 0 -1; -1 -1 0')
Aci = np.linalg.inv(Ac)
Bc = np.matrix('-1 -1 1; 1 0 1; 1 0 -1')
Bci = np.linalg.inv(Bc)
Cc = np.matrix('1 0 1; 0 1 0; 1 0 0')
Cci = np.linalg.inv(Cc)


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
            print(line, end='')
            num, rest = line.split(maxsplit=1)
            num = int(num)
            if rest.startswith('input'):
                ms.append([Layer(num, I, C)])
            elif rest.startswith('gate') or rest.startswith('output'):
                # XXX: watchout!  I'm not sure what'll happen if we have
                # multiple outputs in a circuit
                if rest.startswith('gate'):
                    _, _, arity, _, rest = rest.split(maxsplit=4)
                else:
                    _, _, _, arity, _, rest = rest.split(maxsplit=5)
                arity = int(arity)
                if arity == 1:
                    _, a, b, _, _, _, in1, _ = rest.split(maxsplit=7)
                    a, b, in1 = ints(a, b, in1)
                    if a == 1 and b == 0:
                        # NOT gate
                        a = notgate(ms[in1])
                        ms.append(a)
                    else:
                        raise("error: only support NOT so far:", line.strip())
                elif arity == 2:
                    _, a, b, c, d, _, _, _, in1, in2, _ \
                        = rest.split(maxsplit=10)
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
                        raise("error: only support AND/XOR so far:", line.strip())
                else:
                    raise("error: arity %d unsupported" % arity)
            else:
                raise("error: unknown type")
    return ms[-1]


def eval_bp(bp, inp):
    comp = I
    for m in bp:
        comp = comp * (m.I if inp[m.inp] == '0' else m.J)
    print(comp)
    if (comp == I).all():
        return 0
    elif (comp == C).all():
        return 1
    else:
        print("error: invalid return matrix:\n%s" % comp)
        exit(-1)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("error: invalid arguments")
        exit(-1)

    fname = sys.argv[1]
    inp = sys.argv[2]
    bp = circuit_to_bp(fname)
    print("BRANCHING PROGRAM =", bp)
    out = eval_bp(bp, inp)
    print("OUTPUT = %d" % out)
