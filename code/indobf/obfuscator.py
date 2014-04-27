from __future__ import print_function

from branchingprogram import _group
import _obfuscator as _obf
import utils

from sage.all import copy, VectorSpace, Zmod, ZZ

import copy, os, time
import numpy as np

def to_long(lst):
    return [long(i) for i in lst]

class ObfuscationException(Exception):
    pass

class AbstractObfuscator(object):
    def _print_params(self):
        self.logger('Graded encoding parameters:')
        self.logger('  Security Parameter: %d' % self.secparam)
        self.logger('  Kappa: %d' % self.kappa)
        self.logger('  Alpha: %d' % self.alpha)
        self.logger('  Beta: %d' % self.beta)
        self.logger('  Eta: %d' % self.eta)
        self.logger('  Nu: %d' % self.nu)
        self.logger('  Rho: %d' % self.rho)
        self.logger('  Rho_f: %d' % self.rho_f)
        self.logger('  N: %d' % self.n)

    def _set_params(self, secparam, kappa):
        self.secparam = secparam
        self.kappa = kappa
        self.alpha = self.secparam
        self.beta = self.secparam
        self.rho = self.secparam
        self.rho_f = self.kappa * (self.rho + self.alpha + 2)
        self.eta = self.rho_f + self.alpha + 2 * self.beta + self.secparam + 8
        self.nu = self.eta - self.beta - self.rho_f - self.secparam - 3
        assert self.nu >= self.alpha + self.beta + 5
        if self._use_small_params:
            self.n = self.eta
        else:
            self.n = int(self.eta * np.log2(self.secparam))
        self._print_params()

    def __init__(self, verbose=False, use_small_params=False,
                 use_fast_prime_gen=True):
        self.obfuscation = None
        self._verbose = verbose
        self._use_small_params = use_small_params
        self._use_fast_prime_gen = use_fast_prime_gen
        self.logger = utils.make_logger(self._verbose)
        if self._use_small_params:
            self.logger('* Using small (and *insecure*) parameters for speed')
        if self._use_fast_prime_gen:
            self.logger('* Using CLT13 prime generation')

    def _gen_mlm_params(self, secparam, size, nzs, directory):
        self.logger('Generating MLM parameters...')
        start = time.time()
        if not os.path.exists(directory):
            os.mkdir(directory)
        primes = _obf.setup(secparam, size, self.n, self.alpha, self.beta,
                            self.eta, self.nu, self.rho, nzs, directory,
                            self._use_fast_prime_gen)
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))
        return primes

    def _randomize(self, secparam, bp, primes, islayered):
        self.logger('Randomizing BPs...')
        num = secparam
        start = time.time()
        bps = [copy.deepcopy(bp) for _ in xrange(num)]
        if islayered:
            alphas = None
            for bp, prime in zip(bps, primes[:num]):
                bp.randomize(prime)
        else:
            Rs = [Zmod(prime) for prime in primes[:num]]
            alphas = [[(R.random_element(), R.random_element()) for _ in
                       xrange(len(bp))] for bp, R in zip(bps, Rs)]
            for bp, prime, alpha in zip(bps, primes[:num], alphas):
                bp.randomize(prime, alphas=alpha)
        end = time.time()
        self.logger('Took: %f seconds' % (end - start))
        return bps, alphas

    def _obfuscate(self, bps):
        for i in xrange(len(bps[0])):
            start = time.time()
            self.logger('Obfuscating layer...')
            zeros = [to_long(bp[i].zero.transpose().list()) for bp in bps]
            ones = [to_long(bp[i].one.transpose().list()) for bp in bps]
            _obf.encode_layers(i, bps[0][i].inp, zeros, ones, bps[0][i].zeroset,
                               bps[0][i].oneset)
            end = time.time()
            self.logger('Took: %f seconds' % (end - start))

    def obfuscate(self, bp, secparam, directory):
        start = time.time()
        islayered = True if type(self) == LayeredObfuscator else False
        if bp.randomized:
            raise ObfuscationException('BPs must not be randomized')
        # remove old files in obfuscation directory
        if os.path.isdir(directory):
            files = os.listdir(directory)
            for file in os.listdir(directory):
                path = os.path.join(directory, file)
                os.unlink(path)

        # add two to kappa due to the bookend vectors
        kappa = len(bp) + 2
        # construct straddling sets, and add two to the number of Zs to take
        # bookend vectors into account
        nzs = bp.set_straddling_sets() + 2
        # size is the column/row-length of the matrices
        size = bp.size if islayered else len(_group)

        self._set_params(secparam, kappa)
        self.logger('  Number of Zs: %d' % nzs)
        self.logger('  Size: %d' % size)
        primes = self._gen_mlm_params(secparam, size, nzs, directory)
        bps, alphas = self._randomize(secparam, bp, primes, islayered)
        if not islayered:
            self._construct_multiplicative_constants(bps, alphas)
        self._construct_bookend_vectors(bps, primes, nzs)
        self._obfuscate(bps)
        end = time.time()
        self.logger('Obfuscation took: %f seconds' % (end - start))

    def evaluate(self, directory, inp):
        start = time.time()
        files = os.listdir(directory)
        inputs = sorted(filter(lambda s: 'input' in s, files))
        islayered = True if type(self) == LayeredObfuscator else False
        result = _obf.evaluate(directory, inp, len(inputs), islayered)
        end = time.time()
        self.logger('Evaluating %s took: %f seconds' % (inp, end - start))
        return result

    def encode_benchmark(self):
        _obf.encode_benchmark()

    def cleanup(self):
        _obf.cleanup()

class Obfuscator(AbstractObfuscator):
    def __init__(self, **kwargs):
        super(Obfuscator, self).__init__(**kwargs)

    def _construct_multiplicative_constants(self, bps, alphas):
        start = time.time()
        for layer_idx in xrange(len(bps[0])):
            a0s = []
            a1s = []
            for i in xrange(len(bps)):
                a0, a1 = alphas[i][layer_idx]
                a0s.append([long(a0)])
                a1s.append([long(a1)])
            _obf.encode_scalars(a0s, bps[0][layer_idx].zeroset, "%d.a0_enc" % layer_idx)
            _obf.encode_scalars(a1s, bps[0][layer_idx].oneset, "%d.a1_enc" % layer_idx)
        end = time.time()
        self.logger('Constructing multiplicative constants took: %f seconds' % (end - start))

    def _construct_bookend_vectors(self, bps, primes, nzs):
        start = time.time()
        sidx, tidx = nzs - 2, nzs - 1
        ss = []
        ts = []
        ps = []
        for bp, prime in zip(bps, primes):
            VSZp = VectorSpace(ZZ.residue_field(ZZ.ideal(prime)), _group.length)
            s = VSZp.random_element() * bp.m0i
            t = bp.m0 * VSZp.random_element()
            ss.append(to_long(s))
            ts.append(to_long(t))
            ps.append([long(s * t)])
        _obf.encode_scalars(ps, [sidx, tidx], "p_enc")
        _obf.encode_vectors(ss, [sidx], "s_enc")
        _obf.encode_vectors(ts, [tidx], "t_enc")
        end = time.time()
        self.logger('Constructing bookend vectors took: %f seconds' % (end - start))

class LayeredObfuscator(AbstractObfuscator):
    def __init__(self, **kwargs):
        super(LayeredObfuscator, self).__init__(**kwargs)

    def _construct_bookend_vectors(self, bps, primes, nzs):
        def compute_vectors():
            start = time.time()
            ss = [to_long(bp.e_1 * bp.m0i) for bp in bps]
            ts = [to_long(bp.m0 * bp.e_w) for bp in bps]
            end = time.time()
            self.logger('  Computing bookend vectors: %f seconds' % (end - start))
            return ss, ts
        start = time.time()
        sidx, tidx = nzs - 2, nzs - 1
        ss, ts = compute_vectors()
        _obf.encode_vectors(ss, [sidx], "s_enc")
        _obf.encode_vectors(ts, [tidx], "t_enc")
        end = time.time()
        self.logger('Constructing bookend vectors took: %f seconds' % (end - start))
