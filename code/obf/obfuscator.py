import utils
import os, time

class Obfuscator(object):
    def __init__(self, obf, verbose=False):
        self._state = None
        self._verbose = verbose
        obf.verbose(self._verbose)
        self.logger = utils.make_logger(self._verbose)

    def _remove_old(self, directory):
        # remove old files in obfuscation directory
        if os.path.isdir(directory):
            files = os.listdir(directory)
            for file in os.listdir(directory):
                p = os.path.join(directory, file)
                os.unlink(p)

    def obfuscate(self, circuit, secparam, directory, obliviate=False,
                  nslots=None, kappa=None):
        raise NotImplementedError

    def _evaluate(self, directory, inp, f, obf):
        self.logger('Evaluating %s...' % inp)
        start = time.time()
        files = os.listdir(directory)
        inputs = sorted(filter(lambda s: 'input' in s, files))
        result = f(directory, inp, len(inputs))
        end = time.time()
        self.logger('Took: %f' % (end - start))
        if self._verbose:
            obf.max_mem_usage()
        return result

    def evaluate(self, directory, inp):
        raise NotImplementedError

    def cleanup(self):
        raise NotImplementedError
