from __future__ import print_function

import networkx as nx
from circuit import parse, ParseException
from sage.all import copy, GF, MatrixSpace, VectorSpace, ZZ
import utils

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
    # def conjugate(self, M, Mi):
    #     return Layer(self.inp, Mi * self.zero * M, Mi * self.one * M,
    #                  zeroset=self.zeroset, oneset=self.oneset)
    def group(self, group, prime):
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

class _Graph(object):
    def __init__(self, inp, graph, nlayers, num):
        # function mapping layers to input bits
        self.inp = inp
        # actual graph; each node is a tuple (name, num), where 'name' is the
        # (unique) name of the node, and 'num' is a number
        self.graph = graph
        # total number of layers in the graph
        self.nlayers = nlayers
        self.num = num
    def __len__(self):
        return len(self.graph)
    def __str__(self):
        return repr(self.graph.adj)

def relabel(g, num):
    '''
    Relabels all the nodes in graph 'g' from (name, num') to (name, num).
    '''
    new = [(s, num) for (s, _) in g.nodes()]
    return nx.relabel_nodes(g, dict(zip(g.nodes(), new)))

def relabel_internal(g, num):
    '''
    Relabels only those nodes in graph 'g' which are "internal", which is marked
    by whether they have a '-' in their name.
    '''
    new = []
    for s, _ in g.nodes():
        if '-' in s:
            new.append(('%s-%d' % (s, num), num))
        else:
            new.append((s, num))
    return nx.relabel_nodes(g, dict(zip(g.nodes(), new)))

def relabel_layers(g, min):
    nodes = nx.get_node_attributes(g, 'layer')
    for k, v in nodes.iteritems():
        g.node[k]['layer'] = v + min

def contract(g, a, b, name):
    '''
    Contracts the edge between nodes 'a' and 'b' in graph 'g', replacing them
    with a new node 'name', and returns the new graph.
    '''
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

class BranchingProgram(object):
    def __init__(self, fname, verbose=False, obliviate=False):
        self._verbose = verbose
        self.logger = utils.make_logger(self._verbose)
        self.nlayers = 0
        self.ninputs = None
        self.depth = None
        self.bp = None
        self.randomized = False
        self.zero = None
        self.graph = self._load_formula(fname)
        if obliviate:
            self._obliviate_graph(self.graph)
        self.size = len(self.graph)
        self._to_relaxed_matrix_bp(self.graph)
        if obliviate:
            self.obliviate()
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
    
    def _load_formula(self, fname):
        bp = []
        self.nlayers = 0
        def _new_gate(num):
            #
            # New gates are constructed as follows:
            #          1
            #       ------- accept
            #      /
            # src -
            #      \   0
            #       ------- reject
            #
            # with the source set to Layer 0.
            #
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
            #
            # AND gates are constructed as follows:
            #
            # Given BP_1 and BP_2, we merge the acc node of BP_1 with the src
            # node of BP_2 and the rej node of BP_1 with the rej node of BP_2.
            #
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
            #
            # OR(a, b) = NOT AND(NOT a, NOT b)
            #
            in1not = _not_gate(idx1, bp1, idx1)
            in2not = _not_gate(idx2, bp2, idx2)
            r = _and_gate(num, in1not, idx1, in2not, idx2)
            return _not_gate(num, r, num)
        def _not_gate(num, bp, idx):
            #
            # NOT gates are constructed as follows:
            #
            # Given BP, swap the acc and rej nodes.
            #
            g = nx.relabel_nodes(bp.graph, {('acc', idx): ('rej', idx),
                                            ('rej', idx): ('acc', idx)})
            g = relabel(g, num)
            return _Graph(bp.inp, g, bp.nlayers, num)
        def _xor_gate(num, bp1, idx1, bp2, idx2):
            #
            # XOR gates are constructed as follows:
            #
            # Given BP_1 and BP_2 (where BP_2 is the "smaller" of the two BPs),
            # produce NOT BP_2, merge the acc node of BP_1 with the src node of
            # NOT BP_2, and merge the rej node of BP_1 with the src node of
            # BP_2.
            #
            assert num > idx1 and num > idx2
            if len(bp1.graph) < len(bp2.graph):
                bp1, bp2 = bp2, bp1
            # choose a temporary idx outside the range of possible indices
            tmpidx = len(bp1.graph) + len(bp2.graph)
            t1 = bp1.nlayers
            t2 = bp2.nlayers
            relabel_layers(bp2.graph, t1)
            oldlayer = bp2.graph.node[('src', idx2)]['layer']
            # construct (G_2, not(G_2))
            bp2not = _not_gate(tmpidx, bp2, idx2)
            # Need to relabel the internal wires as bp2not so we don't end up
            # with duplicate node names when we merge the graphs
            bp2not.graph = relabel_internal(bp2not.graph, tmpidx)
            g = nx.union(bp2.graph, bp2not.graph)
            g = contract(g, ('acc', tmpidx), ('acc', idx2), ('acc', num))
            g = contract(g, ('rej', tmpidx), ('rej', idx2), ('rej', num))
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
                    raise Exception("xorgate eval failed on %s!" % inp)
            return _Graph(eval, g, t1 + t2, num)

        gates = {
            'ID': lambda num, in1: _id_gate(num, bp[in1], in1),
            'AND': lambda num, in1, in2: _and_gate(num, bp[in1], in1, bp[in2], in2),
            'OR': lambda num, in1, in2: _or_gate(num, bp[in1], in1, bp[in2], in2),
            'NOT': lambda num, in1: _not_gate(num, bp[in1], in1),
            'XOR': lambda num, in1, in2: _xor_gate(num, bp[in1], in1, bp[in2], in2),
        }
        wires = set()
        def _inp_gate(bp, num):
            bp.append(_new_gate(num))
        def _gate(bp, num, lineno, gate, inputs):
            if wires.intersection(inputs):
                raise ParseException(
                    'Line %d: only Boolean formulas supported' % lineno)
            wires.update(inputs)
            bp.append(gates[gate](num, *inputs))
        bp, self.nlayers, self.ninputs, self.depth = parse(fname, bp, _inp_gate,
                                                           _gate)
        return bp

    def _obliviate_graph(self, graph):
        assert self.ninputs
        # Boolean formulas must have ninputs - 1 gates, and we need to obliviate
        # this such that all gates look like XOR gates (since these gates
        # requires the largest number of vertices).  Thus, we calculate the
        # expected number of vertices and the current number of vertices, and
        # subtract to determine how many dummy vertices we need to add.
        # If we have only one input, there's nothing to obliviate.
        if self.ninputs > 1:
            # One XOR gate requires 5 nodes, and every additional XOR gate adds
            # 2 more nodes
            expected = 5 + 2 * (self.ninputs - 2)
            current = len(graph)
            for i in xrange(expected - current):
                graph.graph.add_node('dummy-%d' % i)
            assert len(graph) == expected

    def _to_relaxed_matrix_bp(self, graph):
        # convert graph to a topological sorting of the vertices
        nodes = nx.topological_sort(graph.graph)
        # the accept vertex must be last
        a = nodes.index(('acc', graph.num))
        b = nodes.index(('rej', graph.num))
        if a < b:
            nodes[b], nodes[a] = nodes[a], nodes[b]
        # the source vertex must be first
        src = nodes.index(('src', graph.num))
        nodes[0], nodes[src] = nodes[src], nodes[0]

        mapping = dict(zip(nodes, range(self.size)))
        g = nx.relabel_nodes(graph.graph, mapping)
        # convert graph to relaxed matrix BP
        self.bp = []
        G = MatrixSpace(GF(2), self.size)
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
        self.zero = G.one()

    def randomize(self, prime):
        assert not self.randomized
        MSZp = MatrixSpace(ZZ.residue_field(ZZ.ideal(prime)), self.size)
        def random_matrix():
            while True:
                m = MSZp.random_element()
                if not m.is_singular() and m.rank() == self.size:
                    return m, m.inverse()
        m0, m0i = random_matrix()
        self.bp[0] = self.bp[0].group(MSZp, prime).mult_left(m0)
        for i in xrange(1, len(self.bp)):
            mi, mii = random_matrix()
            self.bp[i-1] = self.bp[i-1].group(MSZp, prime).mult_right(mii)
            self.bp[i] = self.bp[i].group(MSZp, prime).mult_left(mi)
        self.bp[-1] = self.bp[-1].group(MSZp, prime).mult_right(m0i)
        VSZp = VectorSpace(ZZ.residue_field(ZZ.ideal(prime)), self.size)
        self.s = copy(VSZp.zero())
        self.s[0] = 1
        self.t = copy(VSZp.zero())
        self.t[len(self.t) - 1] = 1
        self.m0, self.m0i = m0, m0i
        self.randomized = True

    def _eval_layered_bp(self, x):
        '''
        Evaluates the layered BP on input bitstring 'x'.
        '''
        assert self.graph
        g = self.graph.graph.copy()
        nodes = nx.get_node_attributes(g, 'layer')
        for layer in xrange(self.nlayers):
            choice = 0 if x[self.graph.inp(layer)] == '0' else 1
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
                    
    def _eval_relaxed_matrix_bp(self, x):
        '''
        Evaluates the relaxed matrix BP on input bitstring 'x'.
        '''
        assert self.bp
        m = self.bp[0]
        comp = m.zero if x[m.inp] == '0' else m.one
        for m in self.bp[1:]:
            comp *= m.zero if x[m.inp] == '0' else m.one
        if self.randomized:
            r = self.s * self.m0i * comp * self.m0 * self.t
            return 1 if r == 1 else 0
        else:
            return 1 if comp[0, comp.nrows() - 1] == 1 else 0

    def evaluate(self, x):
        assert self.bp or self.graph
        if self.bp is None:
            return self._eval_layered_bp(x)
        else:
            return self._eval_relaxed_matrix_bp(x)


