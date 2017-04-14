import pyobf._obfuscator as _obf
from pyobf.obfuscator import Obfuscator
from pyobf.sz_bp import SZBranchingProgram

import os, time

OBFUSCATOR_FLAG_NONE = 0x00
OBFUSCATOR_FLAG_NO_RANDOMIZATION = 0x01
OBFUSCATOR_FLAG_VERBOSE = 0x04

ENCODE_LAYER_RANDOMIZATION_TYPE_NONE = 0x00
ENCODE_LAYER_RANDOMIZATION_TYPE_FIRST = 0x01
ENCODE_LAYER_RANDOMIZATION_TYPE_MIDDLE = 0x02
ENCODE_LAYER_RANDOMIZATION_TYPE_LAST = 0x04

class SZObfuscator(Obfuscator):
    def __init__(self, mmap, base=None, verbose=False, nthreads=None,
                 ncores=None):
        super(SZObfuscator, self).__init__(_obf, mmap, base=base,
                                           verbose=verbose, nthreads=nthreads,
                                           ncores=ncores)

    def _construct_bp(self, fname, formula=True):
        self.logger('Constructing BP...')
        start = time.time()
        bp = SZBranchingProgram(fname, verbose=self._verbose, formula=formula)
        nzs = bp.set_straddling_sets()
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return bp, nzs

    def _init_mmap(self, secparam, kappa, nzs, directory, seed, flags):
        self.logger('Initializing mmap...')
        start = time.time()
        if not os.path.exists(directory):
            os.mkdir(directory)
        self._state = _obf.init(directory, self._mmap, secparam, kappa, nzs,
                                self._nthreads, self._ncores, seed, flags)
        end = time.time()
        self.logger('Took: %f' % (end - start))

    def _obfuscate(self, bp, nzs):
        nencodings = 0
        for i in xrange(len(bp)):
            nrows, ncols = bp[i].matrices[0].shape
            nencodings += nrows * ncols * self._base
        self.logger('Total # Encodings: %d' % nencodings)
        for i in xrange(len(bp)):
            self.logger('Obfuscating layer...')
            mats = [bp[i].matrices[j].tolist() for j in xrange(self._base)]
            nrows, ncols = bp[i].matrices[0].shape
            rflags = ENCODE_LAYER_RANDOMIZATION_TYPE_NONE
            if i == 0:
                rflags |= ENCODE_LAYER_RANDOMIZATION_TYPE_FIRST
            if i == len(bp) - 1:
                rflags |= ENCODE_LAYER_RANDOMIZATION_TYPE_LAST
            if 0 < i < len(bp) - 1:
                rflags |= ENCODE_LAYER_RANDOMIZATION_TYPE_MIDDLE
            pows = [[0] * nzs for _ in xrange(self._base)]
            for j in xrange(self._base):
                for k in bp[i].sets[j]:
                    pows[j][k] = 1
            _obf.encode_layer(self._state, self._base, pows, mats, i, nrows,
                              ncols, bp[i].inp, rflags)

    def obfuscate(self, fname, secparam, directory, kappa=None, formula=True,
                  randomization=True, seed=None):
        start = time.time()
        self._remove_old(directory)
        bp, nzs = self._construct_bp(fname, formula=formula)
        if not kappa:
            kappa = nzs
        flags = OBFUSCATOR_FLAG_NONE
        if self._verbose:
            flags |= OBFUSCATOR_FLAG_VERBOSE
        if not randomization:
            flags |= OBFUSCATOR_FLAG_NO_RANDOMIZATION
        self._init_mmap(secparam, kappa, nzs, directory, seed, flags)
        if self._base is None:
            self._base = len(bp[0].matrices)
        self._obfuscate(bp, nzs)
        _obf.wait(self._state)
        end = time.time()
        self.logger('Obfuscation took: %f' % (end - start))
        if self._verbose:
            _obf.max_mem_usage()

    def evaluate(self, directory, inp):
        base = self._base if self._base is not None else 2
        flags = OBFUSCATOR_FLAG_NONE
        if self._verbose:
            flags |= OBFUSCATOR_FLAG_VERBOSE
        input = [int(i, base) for i in inp]
        result = self._evaluate(directory, input, _obf.evaluate, _obf, flags)
        if self._verbose:
            _obf.max_mem_usage()
        return result

