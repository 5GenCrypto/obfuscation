from __future__ import print_function

import _obfuscator as _obf
import agis_obfuscator, utils
from bp_sz import SZBranchingProgram

import os, time

class SZObfuscator(agis_obfuscator.AGISObfuscator):
    def __init__(self, verbose=False):
        super(SZObfuscator, self).__init__(verbose=verbose)

    def obfuscate(self, circuit, secparam, directory, obliviate=False,
                  nslots=None):
        self._remove_old(directory)
        if nslots is None:
            nslots = secparam

        # create a dummy branching program to determine parameters
        bp = SZBranchingProgram(circuit, verbose=self._verbose,
                                obliviate=obliviate)
        kappa = len(bp)
        nzs = bp.set_straddling_sets()

        start = time.time()
        primes = self._gen_mlm_params(secparam, kappa, 0, nzs, directory)
        bps = self._construct_bps(SZBranchingProgram, nslots, circuit, primes, obliviate)
        self._randomize(secparam, bps, primes)
        self._obfuscate(bps, len(primes))
        end = time.time()
        self.logger('Obfuscation took: %f' % (end - start))
        if self._verbose:
            _obf.max_mem_usage()

    def evaluate(self, directory, inp):
        return self._evaluate(directory, inp, _obf.sz_evaluate)
