#!/usr/bin/env sage -python

from __future__ import print_function
import itertools

from sage.all import *

import groups, utils

def ints(*args):
    return (int(arg) for arg in args)
def flatten(l):
    return list(itertools.chain(*l))

class Layer(object):
    def __init__(self, inp, zero, one):
        self.inp = inp
        self.zero = zero
        self.one = one
        self.zeroset = None
        self.oneset = None
    def __repr__(self):
        return "input: %d\nzero:\n%s\none:\n%s\nzeroset: %s\noneset: %s" % (
            self.inp, self.zero, self.one, self.zeroset, self.oneset)
    def to_raw_string(self):
        return "%d %s %s" % (self.inp, self.zero.numpy().tostring(),
                             self.one.numpy().tostring())
    def conjugate(self, M, Mi):
        return Layer(self.inp, Mi * self.zero * M, Mi * self.one * M)
    def invert(self):
        return Layer(self.inp, self.zero.inverse(), self.one.inverse())
    def group(self, group):
        return Layer(self.inp, group(self.zero), group(self.one))
    def mult_scalar(self, alphas):
        return Layer(self.inp, alphas[0] * self.zero, alphas[1] * self.one)
    def mult_left(self, M):
        return Layer(self.inp, M * self.zero, M * self.one)
    def mult_right(self, M):
        return Layer(self.inp, self.zero * M, self.one * M)

def prepend(layers, M):
    return [layers[0].mult_left(M)] + layers[1:]
def append(layers, M):
    return layers[:-1] + [layers[-1].mult_right(M)]

def conjugate(layers, target, group):
    if len(layers) == 1:
        return [layers[0].conjugate(*group.conjugates[target])]
    else:
        layers = prepend(layers, group.conjugates[target][1])
        return append(layers, group.conjugates[target][0])

def invert(layers):
    if len(layers) == 1:
        return [layers[0].invert()]
    else:
        return [invert(l) for l in reversed(layers)]

def notgate(layers, group):
    if len(layers) == 1:
        return [layers[0].mult_left(group.Cc).mult_right(group.Cc * group.C)]
    else:
        return append(prepend(layers, group.Cc), group.Cc * group.C)

class ParseException(Exception):
    pass

class BranchingProgram(object):
    def __init__(self, fname, type='circuit', verbose=False, group='S5'):
        self._verbose = verbose
        self.logger = utils.make_logger(self._verbose)

        self.bp = None
        self.n_inputs = None
        self.depth = None
        try:
            self.bpgroup = groups.groupmap[group]()
        except KeyError:
            self.logger("Unknown group '%s'.  Defaulting to S5." % group)
            self.bpgroup = groups.S5()
        self._group = self.bpgroup.G
        self.zero = self.bpgroup.I
        self.one = self.bpgroup.C
        self.randomized = False
        self.m0, self.m0i = self._group.one(), self._group.one()

        if type not in ('circuit', 'bp'):
            raise ParseException('invalid type argument')
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
    def __repr__(self):
        return repr(self.bp)

    def parse_param(self, line):
        try:
            _, param, value = line.split()
        except ValueError:
            raise ParseException("Invalid line '%s'" % line)
        param = param.lower()
        try:
            value = int(value)
        except ValueError:
            raise ParseException("Invalid value '%s'" % value)
        if param == 'nins':
            self.n_inputs = value
        elif param == 'depth':
            self.depth = value
        else:
            raise ParseException("Invalid parameter '%s'" % param)

    def load_arity_one_gate(self, rest, bp):
        _, a, b, _, _, _, in1, _ = rest.split(None, 7)
        a, b, in1 = ints(a, b, in1)
        if a == 1 and b == 0:
            # NOT gate
            return notgate(bp[in1], self.bpgroup)
        elif a == 0 and b == 1:
            # identity gate
            return bp[in1]
        else:
            raise ParseException("error: unsupported gate:", rest.strip())

    def load_arity_two_gate(self, rest, bp):
        _, a, b, c, d, _, _, _, in1, in2, _ \
            = rest.split(None, 10)
        a, b, c, d, in1, in2 = ints(a, b, c, d, in1, in2)
        if a == b == c == 0 and d == 1:
            # AND gate
            a = conjugate(bp[in1], 'A', self.bpgroup)
            b = conjugate(bp[in2], 'B', self.bpgroup)
            c = conjugate(bp[in1], 'Ai', self.bpgroup)
            d = conjugate(bp[in2], 'Bi', self.bpgroup)
            return flatten([a, b, c, d])
        elif a == d == 0 and b == c == 1:
            # XOR gate
            if repr(self.bpgroup) == 'S5':
                raise ParseException("XOR gates not supported for S5 group")
            return flatten([bp[in1], bp[in2]])
        else:
            raise ParseException("error: unsupported gate:", rest.strip())

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
                    bp.append([Layer(num, self.bpgroup.I, self.bpgroup.C)])
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
                        raise ParseException('unsupported gate arity %s' % arity)
                else:
                    raise("error: unknown type")
        self.bp = bp[-1]

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
            for i in xrange(self.n_inputs):
                if m.inp == i:
                    newbp.append(m)
                else:
                    newbp.append(Layer(i, self.bpgroup.I, self.bpgroup.I))
        ms_needed = (4 ** self.depth) * self.n_inputs
        for _ in xrange((ms_needed - len(newbp)) // self.n_inputs):
            for i in xrange(self.n_inputs):
                newbp.append(Layer(i, self.bpgroup.I, self.bpgroup.I))
        assert len(newbp) == ms_needed
        self.bp = newbp

    def randomize(self, prime, alphas=None):
        MSZp = MatrixSpace(ZZ.residue_field(ZZ.ideal(prime)), self.bpgroup.length)

        def random_matrix():
            while True:
                m = MSZp.random_element()
                if not m.is_singular():
                    return m, m.inverse()

        m0, m0i = random_matrix()
        self.zero = m0 * MSZp(self.zero) * m0i
        self.one = m0 * MSZp(self.one) * m0i
        self.bp[0] = self.bp[0].group(MSZp).mult_left(m0)
        for i in xrange(1, len(self.bp)):
            mi, mii = random_matrix()
            self.bp[i-1] = self.bp[i-1].group(MSZp).mult_right(mii)
            self.bp[i] = self.bp[i].group(MSZp).mult_left(mi)
        self.bp[-1] = self.bp[-1].group(MSZp).mult_right(m0i)
        self._group = MSZp
        self.m0, self.m0i = m0, m0i
        self.randomized = True

        if alphas:
            assert len(alphas) == len(self.bp)
            for i, alpha in enumerate(alphas):
                self.bp[i] = self.bp[i].mult_scalar(alpha)

    def evaluate(self, inp):
        comp = self._group.identity_matrix()
        for m in self.bp:
            comp = comp * (m.zero if inp[m.inp] == '0' else m.one)
        if comp == self.zero:
            return 0
        elif comp == self.one:
            return 1
        else:
            raise Exception('Evaluation failed!')
