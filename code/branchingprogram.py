#!/usr/bin/env sage -python

import itertools
import sys
from sage.all import *

MATRIX_LENGTH = 3
_P = 12071

# G = SL(MATRIX_LENGTH, GF(3))
G = MatrixSpace(GF(3), MATRIX_LENGTH, MATRIX_LENGTH)
MSZp = MatrixSpace(ZZ.residue_field(ZZ.ideal(_P)), MATRIX_LENGTH, MATRIX_LENGTH)

I = G.one()
A = G([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
B = G([[0, 0, 1], [0, -1, 0], [1, 0, 0]])
C = G([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
Ac = G([[-1, 0, 0], [0, 0, -1], [0, -1, 0]])
Bc = G([[0, 1, 0], [1, 0, 1], [-1, 0, 1]])
Bci = Bc.inverse()
Cc = G([[-1, 0, -1], [0, -1, 0], [0, 0, 1]])

CONJUGATES = {
    'A': (Ac, Ac),
    'B': (Bc, Bci),
    'Ai': (Ac, Ac),
    'Bi': (Bc, Bci)
}

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
    def to_raw_string(self):
        return "%d %s %s" % (self.inp, self.I.numpy().tostring(),
                             self.J.numpy().tostring())
    def conjugate(self, M, Mi):
        return Layer(self.inp, Mi * self.I * M, Mi * self.J * M)
    def invert(self):
        return Layer(self.inp, self.I.inverse(), self.J.inverse())
    def group(self, group):
        return Layer(self.inp, group(self.I), group(self.J))
    def mult_left(self, M):
        return Layer(self.inp, M * self.I, M * self.J)
    def mult_right(self, M):
        return Layer(self.inp, self.I * M, self.J * M)

def prepend(layers, M):
    return [layers[0].mult_left(M)] + layers[1:]
def append(layers, M):
    return layers[:-1] + [layers[-1].mult_right(M)]

def conjugate(layers, target):
    if len(layers) == 1:
        return [layers[0].conjugate(*CONJUGATES[target])]
    else:
        layers = prepend(layers, CONJUGATES[target][1])
        return append(layers, CONJUGATES[target][0])

def invert(layers):
    if len(layers) == 1:
        return [layers[0].invert()]
    else:
        return [invert(l) for l in reversed(layers)]

def notgate(layers):
    if len(layers) == 1:
        return [layers[0].mult_left(Cc).mult_right(Cc * C)]
    else:
        return append(prepend(layers, Cc), Cc * C)

class ParseException(Exception):
    pass

class BranchingProgram(object):
    def __init__(self, fname, type='circuit'):
        self.bp = None
        self.n_inputs = None
        self.depth = None
        self._group = G
        self.I = I
        self.C = C
        if type not in ('circuit', 'bp'):
            raise Exception('invalid type argument')
        if type == 'circuit':
            self.load_circuit(fname)
        else:
            self.load_bp(fname)

    def __len__(self):
        return len(self.bp)

    def __iter__(self):
        return self.bp.__iter__()
    def next(self):
        return self.bp.next()

    def parse_param(self, line):
        try:
            _, param, value = line.split()
        except ValueError:
            raise ParseException("Invalid line '%s'" % line)
        param = param.lower()
        if param not in ('nins', 'depth'):
            raise ParseException("Invalid parameter '%s'" % param)
        try:
            value = int(value)
        except ValueError:
            raise ParseException("Invalid value '%s'" % value)
        if param == 'nins':
            self.n_inputs = value
        elif param == 'depth':
            self.depth = value

    def load_arity_one_gate(self, rest, bp):
        _, a, b, _, _, _, in1, _ = rest.split(None, 7)
        a, b, in1 = ints(a, b, in1)
        if a == 1 and b == 0:
            # NOT gate
            return notgate(bp[in1])
        else:
            raise("error: unsupported gate:", line.strip())

    def load_arity_two_gate(self, rest, bp):
        _, a, b, c, d, _, _, _, in1, in2, _ \
            = rest.split(None, 10)
        a, b, c, d, in1, in2 = ints(a, b, c, d, in1, in2)
        if a == b == c == 0 and d == 1:
            # AND gate
            a = conjugate(bp[in1], 'A')
            b = conjugate(bp[in2], 'B')
            c = conjugate(bp[in1], 'Ai')
            d = conjugate(bp[in2], 'Bi')
            return flatten([a, b, c, d])
        elif a == d == 0 and b == c == 1:
            # XOR gate
            return flatten([bp[in1], bp[in2]])
        else:
            raise("error: unsupported gate:", line.strip())

    def load_circuit(self, fname):
        bp = []
        arity_dict = {
            1: self.load_arity_one_gate,
            2: self.load_arity_two_gate
        }
        with open(fname) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                elif line.startswith(':'):
                    self.parse_param(line)
                    continue
                num, rest = line.split(None, 1)
                num = int(num)
                if rest.startswith('input'):
                    bp.append([Layer(num, I, C)])
                elif rest.startswith('gate') or rest.startswith('output'):
                    # XXX: watchout!  I'm not sure what'll happen if we have
                    # multiple outputs in a circuit
                    if rest.startswith('gate'):
                        _, _, arity, _, rest = rest.split(None, 4)
                    else:
                        _, _, _, arity, _, rest = rest.split(None, 5)
                    try:
                        bp.append(arity_dict[int(arity)](rest, bp))
                    except KeyError:
                        raise Exception('unsupported gate arity %s' % arity)
                else:
                    raise("error: unknown type")
        self.bp = bp[-1]

    def load_bp(self, fname):
        raise NotImplemented()

    def save_bp(self, fname):
        with open(fname, mode='wx') as f:
            if self.n_inputs is not None:
                f.write(': nins %d\n' % self.n_inputs)
            if self.depth is not None:
                f.write(': depth %d\n' % self.depth)
            for m in self.bp:
                f.write("%s\n" % m.to_raw_string())

    def obliviate(self):
        assert self.n_inputs is not None and self.depth is not None
        newbp = []
        for m in self.bp:
            for i in range(self.n_inputs):
                if m.inp == i:
                    newbp.append(m)
                else:
                    newbp.append(Layer(i, I, I))
        ms_needed = (4 ** self.depth) * self.n_inputs
        for _ in range((ms_needed - len(newbp)) // self.n_inputs):
            for i in range(self.n_inputs):
                newbp.append(Layer(i, I, I))
        assert(len(newbp) == ms_needed)
        self.bp = newbp

    def randomize(self):
        def random_matrix():
            while True:
                m = MSZp.random_element()
                if not m.is_singular():
                    return m, m.inverse()
        # TODO: how come multiplying by bookend matrices m0 doesn't work?
        # m0, m0i = random_matrix()
        # self.C = G(m0 * MSZp(C) * m0i)
        self.bp[0] = self.bp[0].group(MSZp)# .mult_left(m0)
        for i in range(1, len(self.bp)):
            mi, mii = random_matrix()
            self.bp[i-1] = self.bp[i-1].mult_right(mii)
            self.bp[i] = self.bp[i].group(MSZp).mult_left(mi)
        self.bp[-1] = self.bp[-1].group(MSZp)# .mult_right(m0i)
        self._group = MSZp

    def evaluate(self, inp):
        comp = self._group.identity_matrix()
        for m in self.bp:
            comp = comp * (m.I if inp[m.inp] == '0' else m.J)
        comp = G(comp)
        if comp == self.I:
            return 0
        elif comp == self.C:
            return 1
        else:
            # raise Exception("invalid return matrix:\n%s" % comp)
            return None
