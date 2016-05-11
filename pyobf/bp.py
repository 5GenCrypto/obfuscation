import pyobf.utils as utils

class Layer(object):
    def __init__(self, inp, matrices, sets):
        self.inp = inp
        self.matrices = matrices
        if sets is None:
            self.sets = [None, None]
        else:
            self.sets = sets
    def size(self):
        size = 0
        for matrix in self.matrices:
            size += len(matrix)
        return size
    def __repr__(self):
        return "\
        input: %d\nzero:\n%s\none:\n%s\nzeroset: %s\noneset: %s" % (
            self.inp, self.matrices[0], self.matrices[1], self.sets[0], self.sets[1])
    def mult_scalar(self, alphas):
        return Layer(self.inp,
                     [alphas[0] * self.matrices[0],
                      alphas[1] * self.matrices[1]],
                     self.sets)
    def mult_left(self, M):
         return Layer(self.inp, [M * self.matrices[0], M * self.matrices[1]],
                      self.sets)
    def mult_right(self, M):
         return Layer(self.inp, [self.matrices[0] * M, self.matrices[1] * M],
                      self.sets)

    
class AbstractBranchingProgram(object):
    def __init__(self, verbose=False):
        self._verbose = verbose
        self.logger = utils.make_logger(self._verbose)
        self.nlayers = 0
        self.ninputs = None
        self.bp = None
        self.zero = None
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
                    layer.sets[0] = [n - 1, n]  if i else [n]
                    layer.sets[1] = [n, n + 1]
                    n += 2
                else:
                    layer.sets[0] = [n - 1, n] if max else [n]
                    layer.sets[1] = [n]
                    n += 1
        return n

    def obliviate(self):
        raise NotImplementedError
    def evaluate(self, x):
        raise NotImplementedError
