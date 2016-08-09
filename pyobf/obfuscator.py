import pyobf.utils as utils
import os, time

MMAP_CLT = 0x00
MMAP_GGHLITE = 0x01
MMAP_DUMMY = 0x02

def get_mmap_flag(mmap):
    if mmap == 'CLT':
        return MMAP_CLT
    elif mmap == 'GGH':
        return MMAP_GGHLITE
    else:
        return MMAP_DUMMY

class Obfuscator(object):
    def __init__(self, obf, mmap, base=None, verbose=False, nthreads=None,
                 ncores=None):
        self._state = None
        self._verbose = verbose
        self._nthreads = nthreads
        self._ncores = ncores
        self._base = base
        self.logger = utils.make_logger(self._verbose)
        self._mmap = get_mmap_flag(mmap)

    def _remove_old(self, directory):
        # remove old files in obfuscation directory
        if os.path.isdir(directory):
            files = os.listdir(directory)
            for file in os.listdir(directory):
                p = os.path.join(directory, file)
                os.unlink(p)

    def obfuscate(self, circuit, secparam, directory, kappa=None, formula=True):
        raise NotImplementedError

    def _evaluate(self, directory, inp, f, obf, flags):
        self.logger('Evaluating %s...' % inp)
        start = time.time()
        files = os.listdir(directory)
        inputs = sorted(filter(lambda s: 'input' in s, files))
        result = f(directory, inp, self._mmap, len(inputs), self._ncores, flags)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        if self._verbose:
            obf.max_mem_usage()
        return result

    def evaluate(self, directory, inp):
        raise NotImplementedError

    def cleanup(self):
        raise NotImplementedError
