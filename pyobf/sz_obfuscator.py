import pyobf._obfuscator as _obf
from pyobf.obfuscator import Obfuscator
from pyobf.sz_bp import SZBranchingProgram

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

    def _gen_mlm_params(self, secparam, kappa, nzs, directory):
        self.logger('Generating MLM parameters...')
        start = time.time()
        if not os.path.exists(directory):
            os.mkdir(directory)
        mlm = 0 if self._mlm == 'CLT' else 1
        self._state = _obf.init(directory, mlm, secparam, kappa, nzs,
                                self._nthreads, self._ncores)
        end = time.time()
        self.logger('Took: %f' % (end - start))

    def _construct_bps(self, fname, obliviate, formula=True):
        self.logger('Constructing BP...')
        start = time.time()
        bp = SZBranchingProgram(fname, verbose=self._verbose, obliviate=obliviate,
                                formula=formula)
        bp.set_straddling_sets()
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return bp

    def _obfuscate(self, bp, nzs):
        for i in xrange(len(bp)):
            self.logger('Obfuscating layer...')
            zeros = bp[i].zero.tolist()
            ones = bp[i].one.tolist()
            nrows, ncols = bp[i].zero.shape
            typ = 0
            if i == 0:
                typ = typ | 1
            if i == len(bp) - 1:
                typ = typ | 4
            if 0 < i < len(bp) - 1:
                typ = typ | 2
            zero_pows = [0] * nzs
            one_pows = [0] * nzs
            for j in bp[i].zeroset:
                zero_pows[j] = 1
            for j in bp[i].oneset:
                one_pows[j] = 1
            _obf.encode_layer(self._state, i, nrows, ncols, bp[i].inp, typ,
                              zero_pows, one_pows, zeros, ones)

    def obfuscate(self, fname, secparam, directory, obliviate=False,
                  kappa=None, formula=True):
        start = time.time()

        self._remove_old(directory)

        # create a dummy branching program to determine parameters
        bp = SZBranchingProgram(fname, verbose=self._verbose,
                                obliviate=obliviate, formula=formula)
        nzs = bp.set_straddling_sets()
        if not kappa:
            kappa = nzs

        self._gen_mlm_params(secparam, kappa, nzs, directory)
        bp = self._construct_bps(fname, obliviate, formula=formula)
        self._obfuscate(bp, nzs)

        _obf.wait(self._state)
        if self._verbose:
            _obf.max_mem_usage()

        end = time.time()
        self.logger('Obfuscation took: %f' % (end - start))

    def evaluate(self, directory, inp):
        return self._evaluate(directory, inp, _obf.evaluate, _obf)
