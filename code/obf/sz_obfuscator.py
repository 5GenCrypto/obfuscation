import _obfuscator as _obf
from agis_obfuscator import AGISObfuscator
from sz_bp import SZBranchingProgram

import time

class SZObfuscator(AGISObfuscator):
    def __init__(self, verbose=False, nthreads=None):
        super(SZObfuscator, self).__init__(verbose=verbose, nthreads=nthreads)

    def obfuscate(self, circuit, secparam, directory, obliviate=False,
                  nslots=None, kappa=None):
        start = time.time()

        self._remove_old(directory)
        if nslots is None:
            nslots = secparam

        # create a dummy branching program to determine parameters
        bp = SZBranchingProgram(circuit, verbose=self._verbose,
                                obliviate=obliviate)
        if not kappa:
            kappa = len(bp)
        nzs = bp.set_straddling_sets()

        primes = self._gen_mlm_params(secparam, kappa, 0, nzs, directory)
        bps = self._construct_bps(SZBranchingProgram, nslots, circuit, primes,
                                  obliviate)
        self._randomize(secparam, bps, primes)
        self._obfuscate(bps, len(primes))

        end = time.time()
        self.logger('Obfuscation took: %f' % (end - start))
        if self._verbose:
            _obf.max_mem_usage()

    def evaluate(self, directory, inp):
        return self._evaluate(directory, inp, _obf.sz_evaluate, _obf)
