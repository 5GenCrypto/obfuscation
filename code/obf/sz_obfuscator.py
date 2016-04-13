import _obfuscator as _obf
from obfuscator import Obfuscator
from sz_bp import SZBranchingProgram

import os, time

def to_long(lst):
    return [long(i) for i in lst]

def pad(array, length, bplength):
    if len(array) < length:
        zeros = [to_long([0] * bplength)]
        return array + (zeros * (length - len(array)))
    else:
        return array

class SZObfuscator(Obfuscator):
    def __init__(self, mlm, verbose=False, nthreads=None, ncores=None):
        super(SZObfuscator, self).__init__(_obf, mlm, verbose=verbose,
                                           nthreads=nthreads, ncores=ncores)

    def _gen_mlm_params(self, secparam, kappa, width, nzs, directory):
        self.logger('Generating MLM parameters...')
        start = time.time()
        if not os.path.exists(directory):
            os.mkdir(directory)
        mlm = 0 if self._mlm == 'CLT' else 1
        self._state, primes = _obf.setup(secparam, kappa, width, nzs, directory,
                                         mlm, self._nthreads, self._ncores)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return primes

    def _construct_bps(self, bpclass, nslots, fname, obliviate, formula=True):
        self.logger('Constructing %d BP...' % nslots)
        start = time.time()
        bps = []
        for _ in xrange(nslots):
            bp = bpclass(fname, verbose=self._verbose, obliviate=obliviate,
                         formula=formula)
            bp.set_straddling_sets()
            bps.append(bp)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return bps

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
            bplength = len(bps[0][i].zero.transpose().list())
            zeros = pad([to_long(bp[i].zero.transpose().list()) for bp in bps],
                        length, bplength)
            ones = pad([to_long(bp[i].one.transpose().list()) for bp in bps],
                       length, bplength)
            nrows = bps[0][i].zero.nrows()
            ncols = bps[0][i].zero.ncols()
            assert(len(bps[0][i].zeroset) == 1)
            assert(len(bps[0][i].oneset) == 1)
            assert(bps[0][i].zeroset[0] == bps[0][i].oneset[0])
            _obf.encode_layers(self._state, i, nrows, ncols, bps[0][i].inp,
                               zeros, ones)

    def obfuscate(self, fname, secparam, directory, obliviate=False,
                  kappa=None, formula=True):
        start = time.time()

        self._remove_old(directory)
        if self._mlm == 'CLT':
            nslots = secparam
        else:
            nslots = 1

        # create a dummy branching program to determine parameters
        bp = SZBranchingProgram(fname, verbose=self._verbose,
                                obliviate=obliviate, formula=formula)
        nzs = bp.set_straddling_sets()
        if not kappa:
            kappa = nzs

        primes = self._gen_mlm_params(secparam, kappa, 0, nzs, directory)
        bps = self._construct_bps(SZBranchingProgram, nslots, fname, obliviate,
                                  formula=formula)
        if self._mlm == 'CLT':
            self._randomize(secparam, bps, primes)
        else:
            print("Warning: No randomization for GGH yet")
        self._obfuscate(bps, len(primes))

        _obf.wait(self._state)
        if self._verbose:
            _obf.max_mem_usage()

        end = time.time()
        self.logger('Obfuscation took: %f' % (end - start))

    def evaluate(self, directory, inp):
        return self._evaluate(directory, inp, _obf.evaluate, _obf)
