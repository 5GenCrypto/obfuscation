from __future__ import print_function

import itertools

from circuit import parse
from sage.all import MatrixSpace, VectorSpace, ZZ
from branchingprogram import AbstractBranchingProgram, ParseException, Layer
import groups

def flatten(l):
    return list(itertools.chain(*l))

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

def notgate(layers, group):
    if len(layers) == 1:
        return [layers[0].mult_left(group.Cc).mult_right(group.Cc * group.C)]
    else:
        return append(prepend(layers, group.Cc), group.Cc * group.C)

class BarringtonBranchingProgram(AbstractBranchingProgram):
    def __init__(self, fname, prime, verbose=False, obliviate=False):
        super(BarringtonBranchingProgram, self).__init__(verbose=verbose)
        self.group = groups.SL3(prime)
        self.zero = self.group.I
        self.one = self.group.C
        self._load_circuit(fname)
        if obliviate:
            self.obliviate()

    def _and_gate(self, in1, in2):
        a = conjugate(in1, 'A', self.group)
        b = conjugate(in2, 'B', self.group)
        c = conjugate(in1, 'Ai', self.group)
        d = conjugate(in2, 'Bi', self.group)
        return flatten([a, b, c, d])

    def _id_gate(self, in1):
        return in1

    def _or_gate(self, in1, in2):
        in1not = self._not_gate(in1)
        in2not = self._not_gate(in2)
        r = self._and_gate(in1not, in2not)
        return self._not_gate(r)

    def _not_gate(self, in1):
        return notgate(in1, self.group)

    def _xor_gate(self, in1, in2):
        if repr(self.group) == 'S5':
            raise ParseException("XOR gates not supported for S5 group")
        return flatten([in1, in2])

    def _load_circuit(self, fname):
        bp = []
        gates = {
            'AND': lambda in1, in2: self._and_gate(bp[in1], bp[in2]),
            'ID': lambda in1: self._id_gate(bp[in1]),
            'OR': lambda in1, in2: self._or_gate(bp[in1], bp[in2]),
            'NOT': lambda in1: self._not_gate(bp[in1]),
            'XOR': lambda in1, in2: self._xor_gate(bp[in1], bp[in2]),
        }
        def _new_gate(bp, num):
            bp.append([Layer(int(num), self.zero, self.one)])
        def _gate(bp, num, lineno, gate, inputs):
            bp.append(gates[gate.upper()](*inputs))
        self.bp, _, self.ninputs, self.depth = parse(fname, bp, _new_gate, _gate)

    def randomize(self, prime, alphas=None):
        assert not self.randomized
        MSZp = MatrixSpace(ZZ.residue_field(ZZ.ideal(prime)), self.group.length)
        def random_matrix():
            while True:
                m = MSZp.random_element()
                if not m.is_singular():
                    return m, m.inverse()
        m0, m0i = random_matrix()
        self.zero = m0 * MSZp(self.zero) * m0i
        self.one = m0 * MSZp(self.one) * m0i
        self.bp[0] = self.bp[0].group(MSZp, prime).mult_left(m0)
        for i in xrange(1, len(self.bp)):
            mi, mii = random_matrix()
            self.bp[i-1] = self.bp[i-1].group(MSZp, prime).mult_right(mii)
            self.bp[i] = self.bp[i].group(MSZp, prime).mult_left(mi)
        self.bp[-1] = self.bp[-1].group(MSZp, prime).mult_right(m0i)
        VSZp = VectorSpace(ZZ.residue_field(ZZ.ideal(prime)), self.group.length)
        self.s = VSZp.random_element()
        self.t = VSZp.random_element()
        self.m0, self.m0i = m0, m0i
        self.randomized = True

        if alphas:
            assert len(alphas) == len(self.bp)
            for i, alpha in enumerate(alphas):
                self.bp[i] = self.bp[i].mult_scalar(alpha)

    def evaluate(self, inp):
        m = self.bp[0]
        comp = m.zero if inp[m.inp] == '0' else m.one
        for m in self.bp[1:]:
            comp *= m.zero if inp[m.inp] == '0' else m.one
        if comp == self.zero:
            return 0
        elif comp == self.one:
            return 1
        else:
            raise Exception('Evaluation failed!')

