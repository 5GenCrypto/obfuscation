from __future__ import print_function

import _obfuscator as _obf
import utils
from bp_sww import SWWBranchingProgram
from bp_barrington import BarringtonBranchingProgram

from sage.all import copy, VectorSpace, Zmod, ZZ
from numpy import log2

import copy, os, time

def to_long(lst):
    return [long(i) for i in lst]

def pad(array, length, bplength):
    if len(array) < length:
        zeros = [to_long([0] * bplength)]
        return array + (zeros * (length - len(array)))
    else:
        return array

class AbstractObfuscator(object):

    def __init__(self, verbose=False):
        self.obfuscation = None
        self._verbose = verbose
        _obf.verbose(self._verbose)
        self.logger = utils.make_logger(self._verbose)

    def _gen_mlm_params(self, secparam, kappa, width, nzs, directory):
        self.logger('Generating MLM parameters...')
        start = time.time()
        if not os.path.exists(directory):
            os.mkdir(directory)
        self._state, primes = _obf.setup(secparam, kappa, width, nzs, directory)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return primes

    def _construct_one_bp(self, circuit, primes, obliviate, bpclass):
        self.logger('Constructing 1 BP (*insecure*)...')
        start = time.time()
        bp = bpclass(circuit, primes[0], verbose=self._verbose,
                     obliviate=obliviate)
        bp.set_straddling_sets()
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return [bp]

    def _construct_bps(self, circuit, primes, obliviate, bpclass):
        self.logger('Constructing %d BPs...' % len(primes))
        start = time.time()
        bps = []
        for prime in primes:
            bp = bpclass(circuit, prime, verbose=self._verbose,
                         obliviate=obliviate)
            bp.set_straddling_sets()
            bps.append(bp)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return bps

    def _randomize(self, secparam, bps, primes, is_sww):
        self.logger('Randomizing BPs...')
        start = time.time()
        if is_sww:
            alphas = None
            for bp, prime in zip(bps, primes):
                bp.randomize(prime)
                bp.set_straddling_sets()
        else:
            rings = [Zmod(prime) for prime in primes]
            alphas = [[(ring.random_element(), ring.random_element()) for _ in
                       xrange(len(bp))] for bp, ring in zip(bps, rings)]
            for bp, prime, alpha in zip(bps, primes, alphas):
                bp.randomize(prime, alpha)
                bp.set_straddling_sets()
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return alphas

    def _obfuscate(self, bps, length):
        for i in xrange(len(bps[0])):
            self.logger('Obfuscating layer...')
            start = time.time()
            bplength = len(bps[0][0].zero.transpose().list())
            zeros = pad([to_long(bp[i].zero.transpose().list()) for bp in bps],
                        length, bplength)
            ones = pad([to_long(bp[i].one.transpose().list()) for bp in bps],
                       length, bplength)
            _obf.encode_layers(self._state, i, bps[0][i].inp, zeros, ones,
                               bps[0][i].zeroset, bps[0][i].oneset)
            end = time.time()
            self.logger('Took: %f' % (end - start))

    def obfuscate(self, circuit, secparam, directory, obliviate=False):
        # remove old files in obfuscation directory
        if os.path.isdir(directory):
            files = os.listdir(directory)
            for file in os.listdir(directory):
                p = os.path.join(directory, file)
                os.unlink(p)

        is_sww = True if type(self) == SWWObfuscator else False
        bpclass = SWWBranchingProgram if is_sww else BarringtonBranchingProgram
        # create a dummy branching program to determine parameters
        bp = bpclass(circuit, 3, verbose=self._verbose, obliviate=obliviate)

        # add two to kappa due to the bookend vectors
        kappa = len(bp) + 2
        # construct straddling sets, and add two to the number of Zs to take
        # bookend vectors into account
        nzs = bp.set_straddling_sets() + 2
        # width is the column/row-length of the matrices
        width = bp.size if is_sww else len(bp.group)

        start = time.time()
        primes = self._gen_mlm_params(secparam, kappa, width, nzs, directory)
        bps = self._construct_one_bp(circuit, primes, obliviate, bpclass)
        # bps = self._construct_bps(circuit, primes, obliviate, bpclass)
        alphas = self._randomize(secparam, bps, primes, is_sww)
        if not is_sww:
            self._construct_multiplicative_constants(bps, alphas)
        self._construct_bookend_vectors(bps, primes, nzs)
        self._obfuscate(bps, len(primes))
        end = time.time()
        self.logger('Obfuscation took: %f' % (end - start))
        if self._verbose:
            _obf.max_mem_usage()

    def evaluate(self, directory, inp):
        self.logger('Evaluating %s...' % inp)
        start = time.time()
        is_sww = True if type(self) == SWWObfuscator else False
        files = os.listdir(directory)
        inputs = sorted(filter(lambda s: 'input' in s, files))
        result = _obf.evaluate(directory, inp, len(inputs), is_sww)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        if self._verbose:
            _obf.max_mem_usage()
        return result

    def attack(self, directory, secparam):
        self.logger('Attacking...')
        start = time.time()
        is_sww = True if type(self) == SWWObfuscator else False
        files = os.listdir(directory)
        inputs = sorted(filter(lambda s: 'input' in s, files))
        bplength = len(inputs)
        kappa = bplength + 2 # add two due to bookend vectors
        nzs = kappa # FIXME:
        inp = '1' * bplength
        result = _obf.attack(directory, inp, len(inputs), is_sww, secparam,
                             kappa, nzs)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return result

    def cleanup(self):
        _obf.cleanup(self._state)


class BarringtonObfuscator(AbstractObfuscator):
    def __init__(self, **kwargs):
        super(BarringtonObfuscator, self).__init__(**kwargs)

    def _construct_multiplicative_constants(self, bps, alphas):
        self.logger('Constructing multiplicative constants...')
        start = time.time()
        for layer_idx in xrange(len(bps[0])):
            a0s = []
            a1s = []
            for i in xrange(len(bps)):
                a0, a1 = alphas[i][layer_idx]
                a0s.append([long(a0)])
                a1s.append([long(a1)])
            _obf.encode_scalars(self._state, a0s, bps[0][layer_idx].zeroset,
                                '%d.a0_enc' % layer_idx)
            _obf.encode_scalars(self._state, a1s, bps[0][layer_idx].oneset,
                                '%d.a1_enc' % layer_idx)
        end = time.time()
        self.logger('Took: %f' % (end - start))

    def _construct_bookend_vectors(self, bps, primes, nzs):
        def compute_vectors():
            start = time.time()
            ss = [to_long(bp.s * bp.m0i) for bp in bps]
            ts = [to_long(bp.m0 * bp.t) for bp in bps]
            ps = [[long(bp.s * bp.t)] for bp in bps]
            end = time.time()
            self.logger('  Computing bookend vectors: %f' % (end - start))
            return ss, ts, ps
        self.logger('Constructing bookend vectors...')
        start = time.time()
        sidx, tidx = nzs - 2, nzs - 1
        ss, ts, ps = compute_vectors()
        _obf.encode_scalars(self._state, ps, [sidx, tidx], 'p_enc')
        _obf.encode_vectors(self._state, ss, [sidx], 's_enc')
        _obf.encode_vectors(self._state, ts, [tidx], 't_enc')
        end = time.time()
        self.logger('Took: %f' % (end - start))


class SWWObfuscator(AbstractObfuscator):
    def __init__(self, **kwargs):
        super(SWWObfuscator, self).__init__(**kwargs)

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
