from __future__ import print_function

from bp import Layer
from circuit import ParseException
from sage.all import matrix, MatrixSpace, ZZ
import utils

import random

def transpose(bps):
    bps.reverse()
    newbps = []
    for bp in bps:
        newbps.append(Layer(bp.inp, bp.zero.transpose(), bp.one.transpose(),
                            bp.zeroset, bp.oneset))
    return newbps

def augment(bps, r):
    def _augment(M, r):
        nrows, ncols = M.dimensions()
        Z_1 = matrix.zero(nrows, r)
        Z_2 = matrix.zero(r, ncols)
        I_r = matrix.identity(r)
        tmp1 = M.augment(Z_1).transpose()
        tmp2 = Z_2.augment(I_r).transpose()
        return tmp1.augment(tmp2).transpose()
    newbps = []
    for bp in bps:
        newbps.append(Layer(bp.inp, _augment(bp.zero, r), _augment(bp.one, r),
                            bp.zeroset, bp.oneset))
    return newbps

def mult_left(bps, m):
    bps[0] = bps[0].mult_left(m)
def mult_right(bps, m):
    bps[-1] = bps[-1].mult_right(m)

class SZBranchingProgram(object):
    def __init__(self, fname, verbose=False, obliviate=False):
        self._verbose = verbose
        self.logger = utils.make_logger(self._verbose)
        self.nlayers = 0
        self.ninputs = None
        self.depth = None
        self.bp = None
        self.randomized = False
        self.zero = None
        self.bp = self._load_formula(fname)
    def __len__(self):
        return len(self.bp)
    def __iter__(self):
        return self.bp.__iter__()
    def next(self):
        return self.bp.next()
    def __getitem__(self, i):
        return self.bp[i]
    def __repr__(self):
        return repr(self.bp)

    def set_straddling_sets(self):
        inpdir = {}
        for layer in self.bp:
            inpdir.setdefault(layer.inp, []).append(layer)
        n = 0
        for layers in inpdir.itervalues():
            max = len(layers) - 1
            for i, layer in enumerate(layers):
                if i < max:
                    layer.zeroset = [n - 1, n]  if i else [n]
                    layer.oneset = [n, n + 1]
                    n += 2
                else:
                    layer.zeroset = [n - 1, n] if max else [n]
                    layer.oneset = [n]
                    n += 1
        return n

    def _load_formula(self, fname):
        def _new_gate(num):
            zero = matrix([1, 0])
            one = matrix([1, 1])
            return [Layer(num, zero, one)]
        def _two_input_gate(bp0, bp1, left, right):
            bp1 = augment(transpose(bp1), 1)
            mult_left(bp1, left)
            mult_right(bp1, right)
            bp0.extend(bp1)
            return bp0
        def _and_gate(num, bp0, bp1):
            left = matrix([[0, 0, 1], [0, 1, 0]])
            right = matrix([[0, 1], [1, 0]])
            return _two_input_gate(bp0, bp1, left, right)
        def _id_gate(num, bp0):
            left = matrix([[0, 0, 1], [1, 0, 0]])
            right = matrix([[0, 1], [1, 0]])
            return _two_input_gate(bp0, bp0, left, right)
        def _or_gate(num, bp0, bp1):
            left = matrix([[0, 1, 1], [1, -1, 0]])
            right = matrix([[0, 1], [1, 0]])
            return _two_input_gate(bp0, bp1, left, right)
        def _not_gate(num, bp0):
            left = matrix([[0, 0, 1], [-1, 0, 0]])
            right = matrix([[0, 1], [1, 1]])
            return _two_input_gate(bp0, bp0, left, right)
        def _xor_gate(num, bp0, bp1):
            left = matrix([[0, 1, 1], [1, -2, 0]])
            right = matrix([[0, 1], [1, 0]])
            return _two_input_gate(bp0, bp1, left, right)
        bp = []
        with open(fname) as f:
            for lineno, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                if line.startswith(':'):
                    continue
                num, rest = line.split(None, 1)
                try:
                    num = int(num)
                except ValueError:
                    raise ParseException(
                        'Line %d: gate index not a number' % lineno)
                gates = {
                    'AND': lambda num, in1, in2: _and_gate(num, bp[in1], bp[in2]),
                    'ID': lambda num, in1: _id_gate(num, bp[in1]),
                    'OR': lambda num, in1, in2: _or_gate(num, bp[in1], bp[in2]),
                    'NOT': lambda num, in1: _not_gate(num, bp[in1]),
                    'XOR': lambda num, in1, in2: _xor_gate(num, bp[in1], bp[in2]),
                }
                if rest.startswith('input'):
                    bp.append(_new_gate(num))
                elif rest.startswith('gate') or rest.startswith('output'):
                    if rest.startswith('output'):
                        output = True
                    _, gate, rest = rest.split(None, 2)
                    inputs = [int(i) for i in rest.split()]
                    try:
                        bp.append(gates[gate](num, *inputs))
                    except KeyError:
                        raise ParseException(
                            'Line %d: unsupported gate %s' % (lineno, gate))
                    except TypeError:
                        raise ParseException(
                            'Line %d: incorrect number of arguments given' % lineno)
        return bp[-1]

    def randomize(self, prime):
        assert not self.randomized
        
        d_0 = self.bp[0].zero.nrows()
        d_1 = self.bp[0].zero.ncols()
        S = matrix.identity(d_0)
        for i in xrange(d_0):
            S[i, i] = random.randint(0, prime - 1)
        MSZp = MatrixSpace(ZZ.residue_field(ZZ.ideal(prime)), d_0, d_1)
        self.bp[0] = self.bp[0].group(MSZp, prime).mult_left(S)
        prev = None
        for i in xrange(2, len(self.bp) - 1):
            d_i = self.bp[i].zero.ncols()
            MSZp = MatrixSpace(ZZ.residue_field(ZZ.ideal(prime)), d_i, d_i)
            cur = MSZp.random_element()
            if prev is not None:
                self.bp[i] = self.bp[i].group(MSZp, prime).mult_left(prev)
            self.bp[i] = self.bp[i].group(MSZp, prime).mult_right(cur)
            prev = cur
        r = self.bp[-1].zero.nrows()
        c = self.bp[-1].zero.ncols()
        T = matrix.identity(c)
        for i in xrange(c):
            T[i, i] = random.randint(0, prime - 1)
        MSZp = MatrixSpace(ZZ.residue_field(ZZ.ideal(prime)), r, c)
        self.bp[-1] = self.bp[-1].group(MSZp, prime).mult_right(T)
        self.randomized = True

    def evaluate(self, x):
        assert self.bp
        m = self.bp[0]
        comp = m.zero if x[m.inp] == '0' else m.one
        for m in self.bp[1:]:
            comp *= m.zero if x[m.inp] == '0' else m.one
        return comp[0, comp.ncols() - 1] != 0
