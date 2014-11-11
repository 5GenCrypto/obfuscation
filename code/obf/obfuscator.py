import utils
import os

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
                  nslots=None):
        raise NotImplementedError

    def evaluate(self, directory, inp):
        raise NotImplementedError

    def cleanup(self):
        raise NotImplementedError
