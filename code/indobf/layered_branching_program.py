from __future__ import print_function

import networkx as nx

from sage.all import copy, GF, MatrixSpace

import utils

class _BranchingProgram(object):
    def __init__(self, inp, graph, nlayers, num):
        self.inp = inp
        self.graph = graph
        self.nlayers = nlayers
        self.num = num

class _Layer(object):
    def __init__(self, inp, zero, one):
        self.inp = inp
        self.zero = zero
        self.one = one
        self.zeroset = None
        self.oneset = None

def relabel(g, num):
    new = [(s, num) for (s, _) in g.nodes()]
    return nx.relabel_nodes(g, dict(zip(g.nodes(), new)))

class ParseException(Exception):
    pass

def contract(g, a, b, name):
    new = 'tmp'
    g.add_node(new)
    for node in g.predecessors(a):
        g.add_edge(node, new, label=g.edge[node][a]['label'])
    for node in g.neighbors(b):
        g.add_edge(new, node, label=g.edge[b][node]['label'])
    for node in g.predecessors(b):
        g.add_edge(node, new, label=g.edge[node][b]['label'])
    g.remove_node(a)
    g.remove_node(b)
    g = nx.relabel_nodes(g, {new: name})
    return g

class LayeredBranchingProgram(object):
    def __init__(self, fname, verbose=False):
        self._verbose = verbose
        self.logger = utils.make_logger(self._verbose)
        self.graph = None
        self.bp = None
        self.nlayers = 0
        self._load_formula(fname)

    def _load_formula(self, fname):
        bp = []
        self.nlayers = 0
        def _new_gate(num):
            g = nx.digraph.DiGraph()
            g.add_node(('src', num), layer=1)
            g.add_node(('acc', num))
            g.add_node(('rej', num))
            g.add_edge(('src', num), ('acc', num), label=1)
            g.add_edge(('src', num), ('rej', num), label=0)
            def eval(inp):
                if inp == 1:
                    return num
                else:
                    raise Exception("eval failed on %s!" % inp)
            return _BranchingProgram(eval, g, 1, num)
        def _and_gate(num, idx1, idx2):
            bp1 = bp[idx1]
            bp2 = bp[idx2]
            t1 = bp1.nlayers
            t2 = bp2.nlayers
            g = nx.union(bp1.graph, bp2.graph)
            g = contract(g, ('acc', idx1), ('src', idx2), ('node-%d' % num, num))
            g = contract(g, ('rej', idx1), ('rej', idx2), ('rej', num))
            g = relabel(g, num)
            g.node[('node-%d' % num, num)]['layer'] = t1 + t2
            def eval(inp):
                if inp <= t1:
                    return bp1.inp(inp)
                elif inp <= t1 + t2:
                    return bp2.inp(t1 + t2 - inp + 1)
                else:
                    raise Exception("eval failed on %s!" % inp)
            return _BranchingProgram(eval, g, t1 + t2, num)
        def _id_gate(num, idx):
            return bp[idx]
        def _not_gate(num, idx):
            bp1 = bp[idx]
            g = nx.relabel_nodes(bp1.graph, {('acc', idx): ('rej', idx),
                                             ('rej', idx): ('acc', idx)})
            g = relabel(g, num)
            return _BranchingProgram(bp1.inp, g, bp1.nlayers, num)
        gates = {
            'ID': _id_gate,
            'AND': _and_gate,
            'NOT': _not_gate,
        }
        output = False
        with open(fname) as f:
            for line in f:
                if line.startswith('#') or line.startswith(':'):
                    continue
                num, rest = line.split(None, 1)
                try:
                    num = int(num)
                except ValueError:
                    raise ParseException("gate index not a number")
                if rest.startswith('input'):
                    bp.append(_new_gate(num))
                    self.nlayers += 1
                elif rest.startswith('gate') or rest.startswith('output'):
                    if rest.startswith('output'):
                        if output:
                            raise ParseException('only support single output gate')
                        else:
                            output = True
                    _, gate, rest = rest.split(None, 2)
                    inputs = [int(i) for i in rest.split()]
                    # try:
                    bp.append(gates[gate.upper()](num, *inputs))
                    # except KeyError:
                    #     raise ParseException("unsupported gate '%s'" % gate)
                    # except TypeError:
                    #     raise ParseException("incorrect number of arguments given")
                else:
                    raise ParseException("unknown type")
        if not output:
            raise ParseException("no output gate found")
        self.graph = bp[-1]
        self.to_relaxed_matrix_bp()

    def to_relaxed_matrix_bp(self):
        g = self.graph.graph
        w = len(g)
        n = self.nlayers
        G = MatrixSpace(GF(2), w)
        nodes = nx.topological_sort(g)
        if nodes.index(('acc', self.graph.num)) != len(nodes) - 1:
            a = nodes.index(('acc', self.graph.num))
            b = nodes.index(('rej', self.graph.num))
            nodes[b], nodes[a] = nodes[a], nodes[b]
        mapping = dict(zip(nodes, range(w)))
        g = nx.relabel_nodes(g, mapping)
        self.bp = []
        for layer in xrange(1, self.nlayers + 1):
            B0 = copy(G.one())
            B1 = copy(G.one())
            for edge in g.edges_iter():
                e = g[edge[0]][edge[1]]
                assert e['label'] in (0, 1)
                if g.node[edge[0]]['layer'] == layer:
                    if e['label'] == 0:
                        B0[edge[0], edge[1]] = 1
                    else:
                        B1[edge[0], edge[1]] = 1
            self.bp.append(_Layer(self.graph.inp, B0, B1))

    def _eval_layered_bp(self, inp):
        assert self.graph is not None
        g = self.graph.graph.copy()
        nodes = nx.get_node_attributes(g, 'layer')
        for layer in xrange(1, self.nlayers + 1):
            choice = 0 if inp[self.graph.inp(layer)] == '0' else 1
            for node in nodes:
                if g.node[node]['layer'] == layer:
                    for neighbor in g.neighbors(node):
                        if g.edge[node][neighbor]['label'] != choice:
                            g.remove_edge(node, neighbor)
        try:
            nx.dijkstra_path(g, ('src', self.graph.num), ('acc', self.graph.num))
            return 1
        except nx.NetworkXNoPath:
            return 0
                    
    def _eval_relaxed_matrix_bp(self, inp):
        assert self.bp is not None
        m = self.bp[0]
        comp = m.zero if inp[m.inp(1)] == '0' else m.one
        for i, m in enumerate(self.bp[1:]):
            comp *= m.zero if inp[m.inp(i + 2)] == '0' else m.one
        if comp[0, comp.nrows() - 1] == 1:
            return 1
        else:
            return 0

    def evaluate(self, inp):
        assert self.bp or self.graph
        if self.bp is None:
            return self._eval_layered_bp(inp)
        else:
            return self._eval_relaxed_matrix_bp(inp)
