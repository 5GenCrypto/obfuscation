from __future__ import print_function

import networkx as nx

import utils

class _BranchingProgram(object):
    def __init__(self, inp, graph, nlayers, num):
        self.inp = inp
        self.graph = graph
        self.nlayers = nlayers
        self.num = num

def relabel(g, num):
    new = [(s, num) for (s, _) in g.nodes()]
    return nx.relabel_nodes(g, dict(zip(g.nodes(), new)))

class ParseException(Exception):
    pass

def contract(g, a, b, num):
    new = (a[0] + '-' + b[0], num)
    g.add_node(new)
    for node in g.predecessors(a):
        g.add_edge(node, new, label=g.edge[node][a]['label'])
    for node in g.neighbors(b):
        g.add_edge(new, node, label=g.edge[b][node]['label'])
    for node in g.predecessors(b):
        g.add_edge(node, new, label=g.edge[node][b]['label'])
    g.remove_node(a)
    g.remove_node(b)
    return g

class LayeredBranchingProgram(object):
    def __init__(self, fname, verbose=False):
        self._verbose = verbose
        self.logger = utils.make_logger(self._verbose)
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
                    raise Exception("eval failed!")
            return _BranchingProgram(eval, g, 1, num)
        def _and_gate(num, idx1, idx2):
            bp1 = bp[idx1]
            bp2 = bp[idx2]
            t1 = bp1.nlayers
            t2 = bp2.nlayers
            g = nx.union(bp1.graph, bp2.graph)
            g = contract(g, ('acc', idx1), ('src', idx2), num)
            g = contract(g, ('rej', idx1), ('rej', idx2), num)
            g = relabel(g, num)
            g.node[('acc-src', num)]['layer'] = self.nlayers
            def eval(inp):
                if inp <= t1:
                    return bp1.inp(inp)
                elif inp <= t1 + t2:
                    return bp2.inp(t1 + t2 - inp + 1)
                else:
                    raise Exception("eval failed!")
            return _BranchingProgram(eval, g, t1 + t2, num)
        def _id_gate(num, idx):
            return bp[idx]
        def _not_gate(num, idx):
            bp1 = bp[idx]
            print(bp1.graph.nodes())
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
        self.bp = bp[-1]

    def evaluate(self, inp):
        # print('Input = %s' % inp)
        g = self.bp.graph.copy()
        nodes = nx.get_node_attributes(g, 'layer')
        for layer in xrange(1, self.nlayers + 1):
            # print('Layer = %s' % layer)
            # print('Input bit = %s' % self.bp.inp(layer))
            choice = 0 if inp[self.bp.inp(layer)] == '0' else 1
            for node in nodes:
                if g.node[node]['layer'] == layer:
                    for neighbor in g.neighbors(node):
                        if g.edge[node][neighbor]['label'] != choice:
                            g.remove_edge(node, neighbor)
        try:
            nx.dijkstra_path(g, ('src', self.bp.num), ('acc', self.bp.num))
            return 1
        except nx.NetworkXNoPath:
            return 0
                    
