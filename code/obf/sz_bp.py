from __future__ import print_function

from bp import AbstractBranchingProgram, Layer
from circuit import ParseException
from sage.all import matrix, MatrixSpace, ZZ

import json, random, sys

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

class SZBranchingProgram(AbstractBranchingProgram):
    def __init__(self, fname, verbose=False, obliviate=False, formula=True):
        super(SZBranchingProgram, self).__init__(verbose=verbose)
        if formula:
            self.bp = self._load_formula(fname)
        else:
            self.bp = self._load_bp(fname)

    def obliviate(self):
        assert self.ninputs and self.depth
        assert not self.randomized
        newbp = []
        for m in self.bp:
            for i in xrange(self.ninputs):
                if m.inp == i:
                    newbp.append(m)
                else:
                    newbp.append(Layer(i, self.zero, self.zero))
        self.bp = newbp

    def _load_bp(self, fname):
        bp = []
        try:
            with open(fname) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    bp_json = json.loads(line)
                    for step in bp_json['steps']:
                        bp.append(
                            Layer(int(step['position']), matrix(step['0']), matrix(step['1'])))

                    assert len(bp_json['outputs'])    == 1 and \
                           len(bp_json['outputs'][0]) == 2
                    first_out = bp_json['outputs'][0][0].lower()
                    if first_out not in ['false', '0']:
                        if first_out not in ['true', '1']:
                            print('warning: interpreting %s as a truthy output' % first_out)
                        bp[-1].zero.swap_columns(0,1)
                        bp[-1].one .swap_columns(0,1)
                    return bp
        except IOError as e:
            print(e)
            sys.exit(1)
        except ValueError as e:
            print('expected numeric position while parsing branching program JSON')
            print(e)
            sys.exit(1)

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
            right = matrix([[1, 0], [0, 1]])
            mult_right(bp0, right)
            return bp0
        def _or_gate(num, bp0, bp1):
            left = matrix([[0, 1, 1], [1, -1, 0]])
            right = matrix([[0, 1], [1, 0]])
            return _two_input_gate(bp0, bp1, left, right)
        def _not_gate(num, bp0):
            right = matrix([[1, 1], [0, -1]])
            mult_right(bp0, right)
            return bp0
        def _xor_gate(num, bp0, bp1):
            left = matrix([[0, 1, 1], [1, -2, 0]])
            right = matrix([[0, 1], [1, 0]])
            return _two_input_gate(bp0, bp1, left, right)
        with open(fname) as f:
            wires = set()
            bp = []
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
                    if wires.intersection(inputs):
                        raise ParseException(
                            'Line %d: only Boolean formulas supported' % lineno)
                    wires.update(inputs)
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
        prev = None
        for i in xrange(0, len(self.bp)):
            print(i)
            nrows = self.bp[i].zero.nrows()
            ncols = self.bp[i].zero.ncols()
            MSZp = MatrixSpace(ZZ.residue_field(ZZ.ideal(prime)), nrows, ncols)
            if i != 0:
                self.bp[i] = self.bp[i].group(MSZp, prime).mult_left(prev.inverse())
            if i != len(self.bp) - 1:
                MSZp_square = MatrixSpace(ZZ.residue_field(ZZ.ideal(prime)), ncols, ncols)
                cur = MSZp_square.random_element()
                self.bp[i] = self.bp[i].group(MSZp, prime).mult_right(cur)
                prev = cur
            print(self.bp[i])
        # compute S * B_0
        # d_0 = self.bp[0].zero.nrows()
        # d_1 = self.bp[0].zero.ncols()
        # S = matrix.identity(d_0)
        # for i in xrange(d_0):
        #     S[i, i] = random.randint(0, prime - 1)
        # MSZp = MatrixSpace(ZZ.residue_field(ZZ.ideal(prime)), d_0, d_1)
        # self.bp[0] = self.bp[0].group(MSZp, prime).mult_left(S)
        # # compute B_ell * T
        # r = self.bp[-1].zero.nrows()
        # c = self.bp[-1].zero.ncols()
        # T = matrix.identity(c)
        # for i in xrange(c):
        #     T[i, i] = random.randint(0, prime - 1)
        # MSZp = MatrixSpace(ZZ.residue_field(ZZ.ideal(prime)), r, c)
        # self.bp[-1] = self.bp[-1].group(MSZp, prime).mult_right(T)
        self.randomized = True

    def evaluate(self, x):
        assert self.bp
        m = self.bp[0]
        comp = m.zero if x[m.inp] == '0' else m.one
        for m in self.bp[1:]:
            comp *= m.zero if x[m.inp] == '0' else m.one
        return comp[0, comp.ncols() - 1] != 0
