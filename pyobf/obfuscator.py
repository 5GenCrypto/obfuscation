import pyobf.utils as utils
import os, time

class Obfuscator(object):
    def __init__(self, obf, mlm, verbose=False, nthreads=None, ncores=None):
        self._state = None
        self._verbose = verbose
        obf.verbose(self._verbose)
        self._nthreads = nthreads
        self._ncores = ncores
        self.logger = utils.make_logger(self._verbose)
        assert mlm in ('CLT', 'GGH')
        self._mlm = mlm

    def _remove_old(self, directory):
        # remove old files in obfuscation directory
        if os.path.isdir(directory):
            files = os.listdir(directory)
            for file in os.listdir(directory):
                p = os.path.join(directory, file)
                os.unlink(p)

    def obfuscate(self, circuit, secparam, directory, obliviate=False,
                  kappa=None, formula=True):
        raise NotImplementedError

    def _evaluate(self, directory, inp, f, obf):
        mlm = 0 if self._mlm == 'CLT' else 1
        self.logger('Evaluating %s...' % inp)
        start = time.time()
        files = os.listdir(directory)
        inputs = sorted(filter(lambda s: 'input' in s, files))
        result = f(directory, inp, mlm, len(inputs), self._ncores)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        if self._verbose:
            obf.max_mem_usage()
        return result

    def evaluate(self, directory, inp):
        raise NotImplementedError

    def cleanup(self):
        raise NotImplementedError
