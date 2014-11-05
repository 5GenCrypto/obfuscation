from __future__ import print_function

import _zobfuscator as _zobf
from circuit import parse
import utils

import networkx as nx
import os, time

class Circuit(object):
    def __init__(self, fname, verbose=False):
        self._verbose = verbose
        self.ninputs = None
        self.depth = None
        self.y_deg = None
        self.x_degs = None
        self.circuit = None
        self._parse(fname)

    def _inp_gate(self, g, num, inp):
        g[0].add_node(num, label=inp)

    def _gate(self, g, num, lineno, gate, inputs):
        g = g[0]
        # def _id_gate(num, x):
        #     g.add_node(num, gate='ID')
        #     g.add_edge(x, num)
        # def _and_gate(num, x, y):
        #     g.add_node(num, gate='AND')
        #     g.add_edge(x, num)
        #     g.add_edge(y, num)
        # def _or_gate(num, x, y):
        #     g.add_node(num, gate='OR')
        #     g.add_edge(x, num)
        #     g.add_edge(y, num)
        # def _not_gate(num, x):
        #     g.add_node(num, gate='NOT')
        #     g.add_edge(x, num)
        # def _xor_gate(num, x, y):
        #     g.add_node(num, gate='XOR')
        #     g.add_edge(x, num)
        #     g.add_edge(y, num)
        def _mul_gate(num, x, y):
            g.add_node(num, gate='MUL')
            g.add_edge(x, num)
            g.add_edge(y, num)
        gates = {
            # 'ID': _id_gate,
            # 'AND': _and_gate,
            # 'OR': _or_gate,
            # 'NOT': _not_gate,
            # 'XOR': _xor_gate,
            'MUL': _mul_gate,
        }
        gates[gate](num, *inputs)
        return [g]

    def _compute_deg(self, circ, inp):
        start = None
        for k, v in circ.node.iteritems():
            if 'label' in v:
                if v['label'] == inp:
                    start = k
        if start is None:
            raise Exception("couldn't find '%s'" % inp)
        path = nx.shortest_path(circ, start, 2) # XXX:
        deg = 0
        for node in path:
            if 'gate' in circ.node[node]:
                if circ.node[node]['gate'] == 'MUL':
                    deg += 1
        print('Degree of %s = %d' % (inp, deg))
        return deg

    def _parse(self, fname):
        g = nx.digraph.DiGraph()
        self.circuit, self.nlayers, \
            self.ninputs, self.depth = parse(fname, [g], self._inp_gate,
                                             self._gate, keyed=True)
        assert self.nlayers and self.ninputs and self.depth
        self.x_degs = [self._compute_deg(self.circuit, 'x')]
        self.y_deg = self._compute_deg(self.circuit, 'y')

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


class ZimmermanObfuscator(object):

    def __init__(self, verbose=False):
        self._state = None
        self._verbose = verbose
        _zobf.verbose(self._verbose)
        self.logger = utils.make_logger(self._verbose)

    def _gen_mlm_params(self, secparam, nzs, directory):
        self.logger('Generating MLM parameters...')
        start = time.time()
        if not os.path.exists(directory):
            os.mkdir(directory)
        self._state = _zobf.setup(secparam, 3, nzs, directory)
        end = time.time()
        self.logger('Took: %f' % (end - start))

    def _obfuscate(self, circ):
        _zobf.encode_circuit(self._state, [0L], 1, 1)

    def obfuscate(self, circuit, secparam, directory, obliviate=False, nslots=None):
        # remove old files in obfuscation directory
        if os.path.isdir(directory):
            files = os.listdir(directory)
            for file in os.listdir(directory):
                p = os.path.join(directory, file)
                os.unlink(p)

        circ = Circuit(circuit)
        nzs = 1 + 4 * sum(circ.x_degs)

        start = time.time()
        self._gen_mlm_params(secparam + circ.depth, nzs, directory)
        self._obfuscate(circ)
        end = time.time()
        self.logger('Obfuscation took: %f' % (end - start))
        # if self._verbose:
        #     _zobf.max_mem_usage()

    def evaluate(self, directory, inp):
        self.logger('Evaluating %s...' % inp)
        start = time.time()
        files = os.listdir(directory)
        inputs = sorted(filter(lambda s: 'input' in s, files))
        result = _zobf.evaluate(directory, inp, 1)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        # if self._verbose:
        #     _zobf.max_mem_usage()
        return result

