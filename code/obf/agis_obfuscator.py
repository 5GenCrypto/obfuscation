from __future__ import print_function

import _obfuscator as _obf
from obfuscator import Obfuscator
from agis_bp import AGISBranchingProgram

import copy, os, time

def to_long(lst):
    return [long(i) for i in lst]

def pad(array, length, bplength):
    if len(array) < length:
        zeros = [to_long([0] * bplength)]
        return array + (zeros * (length - len(array)))
    else:
        return array

class AGISObfuscator(Obfuscator):
    def __init__(self, verbose=False):
        super(AGISObfuscator, self).__init__(_obf, verbose=verbose)

    def _gen_mlm_params(self, secparam, kappa, width, nzs, directory):
        self.logger('Generating MLM parameters...')
        start = time.time()
        if not os.path.exists(directory):
            os.mkdir(directory)
        self._state, primes = _obf.setup(secparam, kappa, width, nzs, directory)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return primes

    def _construct_bps(self, bpclass, nslots, circuit, primes, obliviate):
        self.logger('Constructing %d BP...' % nslots)
        start = time.time()
        bps = []
        for _, prime in zip(xrange(nslots), primes):
            bp = bpclass(circuit, verbose=self._verbose, obliviate=obliviate)
            bp.set_straddling_sets()
            bps.append(bp)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return bps

    def _construct_bookend_vectors(self, bps, primes, nzs):
        def compute_vectors():
            start = time.time()
            length = len(primes)
            bplength = len(bps[0].s)
            ss = pad([to_long(bp.s * bp.m0i) for bp in bps], length, bplength)
            ts = pad([to_long(bp.m0 * bp.t) for bp in bps], length, bplength)
            end = time.time()
            self.logger('  Computing bookend vectors: %f' % (end - start))
            return ss, ts
        self.logger('Constructing bookend vectors...')
        start = time.time()
        sidx, tidx = nzs - 2, nzs - 1
        ss, ts = compute_vectors()
        _obf.encode_vectors(self._state, ss, [sidx], 's_enc')
        _obf.encode_vectors(self._state, ts, [tidx], 't_enc')
        end = time.time()
        self.logger('Took: %f' % (end - start))

    def _randomize(self, secparam, bps, primes):
        self.logger('Randomizing BPs...')
        start = time.time()
        for bp, prime in zip(bps, primes):
            bp.randomize(prime)
            bp.set_straddling_sets()
        end = time.time()
        self.logger('Took: %f' % (end - start))

    def _obfuscate(self, bps, length):
        for i in xrange(len(bps[0])):
            self.logger('Obfuscating layer...')
            start = time.time()
            bplength = len(bps[0][i].zero.transpose().list())
            zeros = pad([to_long(bp[i].zero.transpose().list()) for bp in bps],
                        length, bplength)
            ones = pad([to_long(bp[i].one.transpose().list()) for bp in bps],
                       length, bplength)
            nrows = bps[0][i].zero.nrows()
            ncols = bps[0][i].zero.ncols()
            _obf.encode_layers(self._state, i, nrows, ncols, bps[0][i].inp,
                               zeros, ones, bps[0][i].zeroset, bps[0][i].oneset)
            end = time.time()
            self.logger('Took: %f' % (end - start))

    def obfuscate(self, circuit, secparam, directory, obliviate=False,
                  nslots=None):
        self._remove_old(directory)
        if nslots is None:
            nslots = secparam

        # create a dummy branching program to determine parameters
        bp = AGISBranchingProgram(circuit, verbose=self._verbose,
                                  obliviate=obliviate)

        # add two to kappa due to the bookend vectors
        kappa = len(bp) + 2
        # construct straddling sets, and add two to the number of Zs to take
        # bookend vectors into account
        nzs = bp.set_straddling_sets() + 2
        # width is the column/row-length of the matrices
        width = bp.size

        start = time.time()
        primes = self._gen_mlm_params(secparam, kappa, width, nzs, directory)
        bps = self._construct_bps(AGISBranchingProgram, nslots, circuit,
                                  primes, obliviate)
        self._randomize(secparam, bps, primes)
        self._construct_bookend_vectors(bps, primes, nzs)
        self._obfuscate(bps, len(primes))
        end = time.time()
        self.logger('Obfuscation took: %f' % (end - start))
        if self._verbose:
            _obf.max_mem_usage()

    def _evaluate(self, directory, inp, f):
        self.logger('Evaluating %s...' % inp)
        start = time.time()
        files = os.listdir(directory)
        inputs = sorted(filter(lambda s: 'input' in s, files))
        result = f(directory, inp, len(inputs))
        end = time.time()
        self.logger('Took: %f' % (end - start))
        if self._verbose:
            _obf.max_mem_usage()
        return result

    def evaluate(self, directory, inp):
        return self._evaluate(directory, inp, _obf.evaluate)

    def attack(self, directory, secparam, nslots):
        self.logger('Attacking...')
        start = time.time()
        files = os.listdir(directory)
        inputs = sorted(filter(lambda s: 'input' in s, files))
        bplength = len(inputs)
        kappa = bplength + 2 # add two due to bookend vectors
        result = _obf.attack(directory, len(inputs), secparam, kappa, nslots)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return result

    def cleanup(self):
        _obf.cleanup(self._state)
