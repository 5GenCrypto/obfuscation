from __future__ import print_function
import itertools

from sage.all import MatrixSpace, ZZ

import groups, utils

class ParseException(Exception):
    pass

class AbstractBranchingProgram(object):
    def __init__(self, verbose=False):
        self._verbose = verbose
        self.logger = utils.make_logger(self._verbose)
        self.bp = None
        self.randomized = False
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

    def obliviate(self):
        raise NotImplemented()
    def randomize(self, prime):
        raise NotImplemented()
    def evaluate(self, inp):
        raise NotImplemented()

_group = groups.S6()

def flatten(l):
    return list(itertools.chain(*l))

class Layer(object):
    def __init__(self, inp, zero, one, zeroset=None, oneset=None):
        self.inp = inp
        self.zero = zero
        self.one = one
        self.zeroset = zeroset
        self.oneset = oneset
    def __repr__(self):
        return "input: %d\nzero:\n%s\none:\n%s\nzeroset: %s\noneset: %s" % (
            self.inp, self.zero, self.one, self.zeroset, self.oneset)
    def conjugate(self, M, Mi):
        return Layer(self.inp, Mi * self.zero * M, Mi * self.one * M,
                     zeroset=self.zeroset, oneset=self.oneset)
    def group(self, group):
        return Layer(self.inp, group(self.zero), group(self.one),
                     zeroset=self.zeroset, oneset=self.oneset)
    def mult_scalar(self, alphas):
        return Layer(self.inp, alphas[0] * self.zero, alphas[1] * self.one,
                     zeroset=self.zeroset, oneset=self.oneset)
    def mult_left(self, M):
        return Layer(self.inp, M * self.zero, M * self.one,
                     zeroset=self.zeroset, oneset=self.oneset)
    def mult_right(self, M):
        return Layer(self.inp, self.zero * M, self.one * M,
                     zeroset=self.zeroset, oneset=self.oneset)

def prepend(layers, M):
    return [layers[0].mult_left(M)] + layers[1:]
def append(layers, M):
    return layers[:-1] + [layers[-1].mult_right(M)]

def conjugate(layers, target):
    if len(layers) == 1:
        return [layers[0].conjugate(*_group.conjugates[target])]
    else:
        layers = prepend(layers, _group.conjugates[target][1])
        return append(layers, _group.conjugates[target][0])

def notgate(layers):
    if len(layers) == 1:
        return [layers[0].mult_left(_group.Cc).mult_right(_group.Cc * _group.C)]
    else:
        return append(prepend(layers, _group.Cc), _group.Cc * _group.C)

class BranchingProgram(AbstractBranchingProgram):
    def __init__(self, fname, verbose=False):
        super(BranchingProgram, self).__init__(verbose=verbose)
        self.n_inputs = None
        self.depth = None
        self.zero = _group.I
        self.one = _group.C
        self._load_circuit(fname)

    def _and_gate(self, in1, in2):
        a = conjugate(in1, 'A')
        b = conjugate(in2, 'B')
        c = conjugate(in1, 'Ai')
        d = conjugate(in2, 'Bi')
        return flatten([a, b, c, d])

    def _id_gate(self, in1):
        return in1

    def _or_gate(self, in1, in2):
        in1not = self._not_gate(in1)
        in2not = self._not_gate(in2)
        r = self._and_gate(in1not, in2not)
        return self._not_gate(r)

    def _not_gate(self, in1):
        return notgate(in1)

    def _xor_gate(self, in1, in2):
        if repr(_group) == 'S5':
            raise ParseException("XOR gates not supported for S5 group")
        return flatten([in1, in2])

    def _parse_param(self, line):
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

    def _load_circuit(self, fname):
        bp = []
        gates = {
            'AND': lambda in1, in2: self._and_gate(bp[in1], bp[in2]),
            'ID': lambda in1: self._id_gate(bp[in1]),
            'OR': lambda in1, in2: self._or_gate(bp[in1], bp[in2]),
            'NOT': lambda in1: self._not_gate(bp[in1]),
            'XOR': lambda in1, in2: self._xor_gate(bp[in1], bp[in2]),
        }
        output = False
        with open(fname) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                elif line.startswith(':'):
                    self._parse_param(line)
                    continue
                num, rest = line.split(None, 1)
                if rest.startswith('input'):
                    bp.append([Layer(int(num), self.zero, self.one)])
                elif rest.startswith('gate') or rest.startswith('output'):
                    if rest.startswith('output'):
                        if output:
                            raise ParseException("only support single output gate")
                        else:
                            output = True
                    _, gate, rest = rest.split(None, 2)
                    inputs = [int(i) for i in rest.split()]
                    try:
                        bp.append(gates[gate.upper()](*inputs))
                    except KeyError:
                        raise ParseException("unsupported gate '%s'" % gate)
                    except TypeError:
                        raise ParseException("incorrect number of arguments given")
                else:
                    raise ParseException("unknown type")
        if not output:
            raise ParseException("no output gate found")
        self.bp = bp[-1]

    def obliviate(self):
        assert self.n_inputs and self.depth
        assert not self.randomized
        newbp = []
        for m in self.bp:
            for i in xrange(self.n_inputs):
                if m.inp == i:
                    newbp.append(m)
                else:
                    newbp.append(Layer(i, self.zero, self.zero))
        ms_needed = (4 ** self.depth) * self.n_inputs
        for _ in xrange((ms_needed - len(newbp)) // self.n_inputs):
            for i in xrange(self.n_inputs):
                newbp.append(Layer(i, self.zero, self.zero))
        assert len(newbp) == ms_needed
        self.bp = newbp

    def randomize(self, prime, alphas=None):
        assert not self.randomized

        MSZp = MatrixSpace(ZZ.residue_field(ZZ.ideal(prime)), _group.length)

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
