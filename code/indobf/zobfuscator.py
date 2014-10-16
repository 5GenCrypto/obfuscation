from __future__ import print_function

from circuit import parse
import networkx as nx
from uuid import uuid4 as uuid

class Circuit(object):
    def __init__(self, fname, verbose=False):
        self._verbose = verbose
        self.ninputs = None
        self.depth = None
        self.y_deg = None
        self.x_degs = None
        self.circuit = None
        self._parse(fname)

    def _inp_gate(self, g, num):
        g[0].add_node(num)

    def _gate(self, g, num, lineno, gate, inputs):
        g = g[0]
        def _id_gate(num, x):
            g.add_node(num, gate='ID')
            g.add_edge(x, num)
        def _and_gate(num, x, y):
            g.add_node(num, gate='AND')
            g.add_edge(x, num)
            g.add_edge(y, num)
        def _or_gate(num, x, y):
            g.add_node(num, gate='OR')
            g.add_edge(x, num)
            g.add_edge(y, num)
        def _not_gate(num, x):
            g.add_node(num, gate='NOT')
            g.add_edge(x, num)
        def _xor_gate(num, x, y):
            g.add_node(num, gate='XOR')
            g.add_edge(x, num)
            g.add_edge(y, num)
        gates = {
            'ID': _id_gate,
            'AND': _and_gate,
            'OR': _or_gate,
            'NOT': _not_gate,
            'XOR': _xor_gate,
        }
        gates[gate](num, *inputs)
        return [g]

    def _parse(self, fname):
        g = nx.digraph.DiGraph()
        self.circuit, self.nlayers, \
            self.ninputs, self.depth = parse(fname, [g], self._inp_gate, self._gate)
        assert self.nlayers and self.ninputs and self.depth
        # TODO:
        # - compute deg(y)
        # - compute deg(x_i) for all i

    def evaluate(self, x):
        assert self.circuit
        g = self.circuit.copy()
        for node in nx.topological_sort(g):
            if node < self.ninputs:
                g.add_node(node, value=int(x[node]))
            elif g.node[node]['gate'] in ('ID', 'NOT'):
                idx = g.pred[node].keys()[0]
                if g.node[node]['gate'] == 'ID':
                    value = g.node[idx]['value']
                elif g.node[node]['gate'] == 'NOT':
                    value = int(g.node[idx]['value'] == 0)
                g.add_node(node, value=value)
            elif g.node[node]['gate'] in ('AND', 'OR', 'XOR'):
                idx1 = g.pred[node].keys()[0]
                idx2 = g.pred[node].keys()[1]
                if g.node[node]['gate'] == 'AND':
                    value = g.node[idx1]['value'] & g.node[idx2]['value']
                elif g.node[node]['gate'] == 'OR':
                    value = g.node[idx1]['value'] | g.node[idx2]['value']
                elif g.node[node]['gate'] == 'XOR':
                    value = g.node[idx1]['value'] ^ g.node[idx2]['value']
                g.add_node(node, value=value)
            else:
                raise Exception('Unable to evaluate')
        idx = nx.topological_sort(g)[-1]
        return g.node[idx]['value']


def straddling_sets(n):
    '''
    Creates an 'n'-straddling set system.
    '''
    symbols = [str(uuid()) for _ in xrange(2 * n - 1)]
    set1, set2 = [], []
    tmp1, tmp2 = set(), set()
    for i, symbol in enumerate(symbols, 1):
        tmp1.add(symbol)
        tmp2.add(symbol)
        if (i - 1) % 2 == 0: # i is odd
            set1.append(tmp1)
            tmp1 = set()
        else:
            set2.append(tmp2)
            tmp2 = set()
    set2.append(tmp2)
    return set1, set2

class ZimmermanObfuscator(object):

    def __init__(self, verbose=False):
        self._state = None
        self._verbose = verbose
        _zobf.verbose(self._verbose)
        self.logger = utils.make_logger(self._verbose)

    def _gen_mlm_params(self, secparam, size, nzs, directory):
        self.logger('Generating MLM parameters...')
        start = time.time()
        if not os.path.exists(directory):
            os.mkdir(directory)
        self._state, primes = _obf.setup(secparam, size, 0, nzs, directory)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return primes

    def obfuscate(self, circuit, secparam, directory):
        # remove old files in obfuscation directory
        if os.path.isdir(directory):
            files = os.listdir(directory)
            for file in os.listdir(directory):
                p = os.path.join(directory, file)
                os.unlink(p)

        circ = Circuit(circuit)
        # set0, set1 = straddling_sets(circ.ninputs)
        nzs = circ.y_deg + 2 * sum(circ.x_degs) + 3 * circ.ninputs

        start = time.time()
        primes = self._gen_mlm_params(secparam + circ.depth, nzs, directory)
        self._obfuscate(circ, primes)
        end = time.time()
        self.logger('Obfuscation took: %f' % (end - start))
        if self._verbose:
            _obf.max_mem_usage()

