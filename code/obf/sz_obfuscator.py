import _obfuscator as _obf
from agis_obfuscator import AGISObfuscator
from sz_bp import SZBranchingProgram

import time

class SZObfuscator(AGISObfuscator):
    def __init__(self, mlm, verbose=False, nthreads=None, ncores=None):
        super(SZObfuscator, self).__init__(mlm, verbose=verbose,
                                           nthreads=nthreads, ncores=ncores)

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
        return self._evaluate(directory, inp, _obf.sz_evaluate, _obf)
