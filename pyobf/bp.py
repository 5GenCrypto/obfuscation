import pyobf.utils as utils

class Layer(object):
    def __init__(self, inp, matrices, sets):
        self.inp = inp
        self.matrices = matrices
        if sets is None:
            self.sets = [None] * len(matrices)
        else:
            self.sets = sets
    def size(self):
        size = 0
        for matrix in self.matrices:
            size += len(matrix)
        return size
    def __repr__(self):
        str = "\ninput: %d\n" % self.inp
        for i, mat in enumerate(self.matrices):
            str += "%d-mat:\n%s\n" % (i, mat)
        for i, set in enumerate(self.sets):
            str += "%d-set: %s\n" % (i, set)
        return str
    def mult_scalar(self, alphas):
        mats = [alphas[i] * self.matrices[i] for i in len(self.matrices)]
        return Layer(self.inp, mats, self.sets)
    def mult_left(self, M):
        mats = [M * mat for mat in self.matrices]
        return Layer(self.inp, mats, self.sets)
    def mult_right(self, M):
        mats = [mat * M for mat in self.matrices]
        return Layer(self.inp, mats, self.sets)


class AbstractBranchingProgram(object):
    def __init__(self, base=None, verbose=False):
        self._verbose = verbose
        self.logger = utils.make_logger(self._verbose)
        self.nlayers = 0
        self.ninputs = None
        self.bp = None
        self.base = base
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
            if len(layers) == 1:
                for i in xrange(len(layers[0].sets)):
                    layers[0].sets[i] = [n]
                n += 1
            else:
                if self.base is not None and self.base != 2:
                    print("Error: straddling sets do not work for base != 2")
                    raise NotImplementedError
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

    def evaluate(self, x):
        raise NotImplementedError
