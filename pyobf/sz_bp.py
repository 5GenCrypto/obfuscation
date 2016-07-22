from __future__ import print_function

from pyobf.bp import AbstractBranchingProgram, Layer
from pyobf.circuit import ParseException

import numpy as np
from numpy import matrix

import json, random, sys

def transpose(bps):
    bps.reverse()
    newbps = []
    for bp in bps:
        newbps.append(Layer(bp.inp,
                            [bp.matrices[0].transpose(),
                             bp.matrices[1].transpose()],
                            bp.sets))
    return newbps

def augment(bps, r):
    def _augment(M, r):
        nrows, ncols = M.shape
        Z_1 = np.zeros([nrows, r], int)
        Z_2 = np.zeros([r, ncols], int)
        I_r = np.identity(r, int)
        tmp1 = np.concatenate((M, Z_1), 1).transpose()
        tmp2 = np.concatenate((Z_2, I_r), 1).transpose()
        return np.concatenate((tmp1, tmp2), 1).transpose()
    newbps = []
    for bp in bps:
        newbps.append(Layer(bp.inp,
                            [_augment(bp.matrices[0], r),
                             _augment(bp.matrices[1], r)],
                            bp.sets))
    return newbps

def mult_left(bps, m):
    bps[0] = bps[0].mult_left(m)
def mult_right(bps, m):
    bps[-1] = bps[-1].mult_right(m)

def swap_columns(m, a, b):
    col = m[:,a].copy()
    m[:,a] = m[:,b].copy()
    m[:,b] = col

class SZBranchingProgram(AbstractBranchingProgram):
    def __init__(self, fname, verbose=False, obliviate=False, formula=True):
        super(SZBranchingProgram, self).__init__(verbose=verbose)
        if formula:
            self.bp = self._load_formula(fname)
        else:
            self.bp = self._load_bp(fname)

    def _load_bp(self, fname):
        bp = []
        try:
            with open(fname) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    bp_json = json.loads(line)
                    for step in bp_json['steps']:
                        keys = list(step)
                        keys = filter(lambda x: x != 'position', keys)
                        keys = sorted(keys)
                        bp.append(
                            Layer(int(step['position']),
                                  [matrix(step[key]) for key in keys],
                                  None))

                    assert len(bp_json['outputs'])    == 1 and \
                           len(bp_json['outputs'][0]) == 2
                    first_out = bp_json['outputs'][0][0].lower()
                    if first_out not in ['false', '0']:
                        if first_out not in ['true', '1']:
                            print('warning: interpreting %s as a truthy output' % first_out)
                        for i in range(len(keys)):
                            swap_columns(bp[-1].matrices[i], 0, 1)
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
            return [Layer(num, [zero, one], None)]
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
        def _parse_file(f):
            wires = set()
            bp = []
            for lineno, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                if line.startswith(':'):
                    continue
                try:
                    num, rest = line.split(None, 1)
                except ValueError:
                    raise ParseException(
                        'Line %d: unable to parse line' % lineno)
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
        try:
            with open(fname) as f:
                return _parse_file(f)
        except IOError as err:
            raise ParseException(err)

    def evaluate(self, x):
        assert self.bp
        m = self.bp[0]
        comp = m.matrices[0] if x[m.inp] == '0' else m.matrices[1]
        for m in self.bp[1:]:
            comp = np.dot(comp, m.matrices[0] if x[m.inp] == '0' else m.matrices[1])
        return comp[0, comp.shape[1] - 1] != 0
