from __future__ import print_function

import networkx as nx

from sage.all import copy, GF, MatrixSpace, VectorSpace, ZZ

from branchingprogram import AbstractBranchingProgram, ParseException, Layer
import utils

class _Graph(object):
    def __init__(self, inp, graph, nlayers, num):
        self.inp = inp
        self.graph = graph
        self.nlayers = nlayers
        self.num = num

def relabel(g, num):
    new = [(s, num) for (s, _) in g.nodes()]
    return nx.relabel_nodes(g, dict(zip(g.nodes(), new)))

def relabel_layers(g, min):
    nodes = nx.get_node_attributes(g, 'layer')
    for k, v in nodes.iteritems():
        g.node[k]['layer'] = v + min

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

class LayeredBranchingProgram(AbstractBranchingProgram):
    def __init__(self, fname, verbose=False):
        super(LayeredBranchingProgram, self).__init__(verbose=verbose)
        self.graph = None
        self.size = None
        self.zero = None
        self.nlayers = 0
        self._load_formula(fname)

    def _load_formula(self, fname):
        # XXX: wow, this code is a total hack!  good luck trying to understand
        # it...
        bp = []
        self.nlayers = 0
        def _new_gate(num):
            g = nx.digraph.DiGraph()
            g.add_node(('src', num), layer=0)
            g.add_node(('acc', num))
            g.add_node(('rej', num))
            g.add_edge(('src', num), ('acc', num), label=1)
            g.add_edge(('src', num), ('rej', num), label=0)
            def eval(inp):
                if inp == 0:
                    return num
                else:
                    raise Exception("newgate eval failed on %s!" % inp)
            return _Graph(eval, g, 1, num)
        def _and_gate(num, bp1, idx1, bp2, idx2):
            t1 = bp1.nlayers
            t2 = bp2.nlayers
            relabel_layers(bp2.graph, t1)
            oldlayer = bp2.graph.node[('src', idx2)]['layer']
            newnode = ('node-%d' % num, num)
            g = nx.union(bp1.graph, bp2.graph)
            g = contract(g, ('acc', idx1), ('src', idx2), newnode)
            g = contract(g, ('rej', idx1), ('rej', idx2), ('rej', num))
            g = relabel(g, num)
            g.node[newnode]['layer'] = oldlayer
            def eval(inp):
                if inp <= t1 - 1:
                    return bp1.inp(inp)
                elif inp <= t1 + t2 - 1:
                    return bp2.inp(inp - t1)
                else:
                    raise Exception("andgate eval failed on %s!" % inp)
            return _Graph(eval, g, t1 + t2, num)
        def _id_gate(num, bp, idx):
            bp.graph = relabel(bp.graph, num)
            bp.num = num
            return bp
        def _or_gate(num, bp1, idx1, bp2, idx2):
            in1not = _not_gate(idx1, bp1, idx1)
            in2not = _not_gate(idx2, bp2, idx2)
            r = _and_gate(num, in1not, idx1, in2not, idx2)
            return _not_gate(num, r, num)
        def _not_gate(num, bp, idx):
            g = nx.relabel_nodes(bp.graph, {('acc', idx): ('rej', idx),
                                            ('rej', idx): ('acc', idx)})
            g = relabel(g, num)
            return _Graph(bp.inp, g, bp.nlayers, num)
        def _xor_gate(num, bp1, idx1, bp2, idx2):
            assert num > idx1 and num > idx2
            assert len(bp1.graph) >= len(bp2.graph)
            tmpidx = len(bp1.graph) + len(bp2.graph)
            t1 = bp1.nlayers
            t2 = bp2.nlayers
            relabel_layers(bp2.graph, t1)
            oldlayer = bp2.graph.node[('src', idx2)]['layer']
            # construct (G_2, not(G_2))
            bp2not = _not_gate(tmpidx, bp2, idx2)
            g = nx.union(bp2.graph, bp2not.graph)
            g.add_edge(('acc', idx2), ('acc', idx1), label=0)
            g.add_edge(('rej', idx1), ('rej', idx2), label=0)
            g = contract(g, ('acc', idx2), ('acc', idx1), ('acc', num))
            g = contract(g, ('rej', idx1), ('rej', idx2), ('rej', num))
            # construct XOR(G_1, G_2)
            g = nx.union(bp1.graph, g)
            accnode = ('acc-%d' % num, num)
            rejnode = ('rej-%d' % num, num)
            g = contract(g, ('acc', idx1), ('src', tmpidx), accnode)
            g = contract(g, ('rej', idx1), ('src', idx2), rejnode)
            g = relabel(g, num)
            g.node[accnode]['layer'] = oldlayer
            g.node[rejnode]['layer'] = oldlayer
            def eval(inp):
                if inp <= t1 - 1:
                    return bp1.inp(inp)
                elif inp <= t1 + t2 - 1:
                    return bp2.inp(inp - t1)
                else:
                    raise Exception("andgate eval failed on %s!" % inp)
            return _Graph(eval, g, t1 + t2, num)

        gates = {
            'ID': lambda num, in1: _id_gate(num, bp[in1], in1),
            'AND': lambda num, in1, in2: _and_gate(num, bp[in1], in1, bp[in2], in2),
            'OR': lambda num, in1, in2: _or_gate(num, bp[in1], in1, bp[in2], in2),
            'NOT': lambda num, in1: _not_gate(num, bp[in1], in1),
            'XOR': lambda num, in1, in2: _xor_gate(num, bp[in1], in1, bp[in2], in2),
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
                    try:
                        bp.append(gates[gate.upper()](num, *inputs))
                    except KeyError:
                        raise ParseException("unsupported gate '%s'" % gate)
                    except TypeError:
                        raise ParseException("incorrect number of arguments given")
                else:
                    raise ParseException("unknown type")
        if not output:
            raise ParseException("no output gate found")
        self.graph = bp[-1]
        # self._obliviate_graph(self.graph)
        self.size = len(self.graph.graph)
        self._to_relaxed_matrix_bp(self.graph)

    def _obliviate_graph(self, graph):
        graph.graph.add_node(('dummy', 0))

    def _to_relaxed_matrix_bp(self, graph):
        # convert graph to a topological sorting of the vertices, with the
        # accept node coming last
        nodes = nx.topological_sort(graph.graph)
        # if nodes.index(('acc', graph.num)) != len(nodes) - 1:
        #     a = nodes.index(('acc', graph.num))
        #     b = nodes.index(('rej', graph.num))
        #     nodes[b], nodes[a] = nodes[a], nodes[b]
        a = nodes.index(('acc', graph.num))
        b = nodes.index(('rej', graph.num))
        if a < b:
            nodes[b], nodes[a] = nodes[a], nodes[b]
        mapping = dict(zip(nodes, range(self.size)))
        g = nx.relabel_nodes(graph.graph, mapping)
        # convert graph to relaxed matrix BP
        self.bp = []
        G = MatrixSpace(GF(2), self.size)
        self.zero = G.one()
        for layer in xrange(self.nlayers):
            zero = copy(G.one())
            one = copy(G.one())
            for edge in g.edges_iter():
                e = g[edge[0]][edge[1]]
                assert e['label'] in (0, 1)
                if g.node[edge[0]]['layer'] == layer:
                    if e['label'] == 0:
                        zero[edge[0], edge[1]] = 1
                    else:
                        one[edge[0], edge[1]] = 1
            self.bp.append(Layer(graph.inp(layer), zero, one))

    def randomize(self, prime):
        assert not self.randomized
        MSZp = MatrixSpace(ZZ.residue_field(ZZ.ideal(prime)), self.size)
        def random_matrix():
            while True:
                m = MSZp.random_element()
                if not m.is_singular() and m.rank() == self.size:
                    return m, m.inverse()
        m0, m0i = random_matrix()
        self.bp[0] = self.bp[0].group(MSZp).mult_left(m0)
        for i in xrange(1, len(self.bp)):
            mi, mii = random_matrix()
            self.bp[i-1] = self.bp[i-1].group(MSZp).mult_right(mii)
            self.bp[i] = self.bp[i].group(MSZp).mult_left(mi)
        self.bp[-1] = self.bp[-1].group(MSZp).mult_right(m0i)
        VSZp = VectorSpace(ZZ.residue_field(ZZ.ideal(prime)), self.size)
        self.e_1 = copy(VSZp.zero())
        self.e_1[0] = 1
        self.e_w = copy(VSZp.zero())
        self.e_w[len(self.e_w) - 1] = 1
        self.m0, self.m0i = m0, m0i
        self.randomized = True

    def _eval_layered_bp(self, inp):
        assert self.graph
        g = self.graph.graph.copy()
        nodes = nx.get_node_attributes(g, 'layer')
        for layer in xrange( self.nlayers):
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
        assert self.bp
        m = self.bp[0]
        comp = m.zero if inp[m.inp] == '0' else m.one
        for i, m in enumerate(self.bp[1:]):
            comp *= m.zero if inp[m.inp] == '0' else m.one
        if self.randomized:
            r = self.e_1 * self.m0i * comp * self.m0 * self.e_w
            return 1 if r == 1 else 0
        else:
            return 1 if comp[0, comp.nrows() - 1] == 1 else 0

    def evaluate(self, inp):
        assert self.bp or self.graph
        if self.bp is None:
            return self._eval_layered_bp(inp)
        else:
            return self._eval_relaxed_matrix_bp(inp)
