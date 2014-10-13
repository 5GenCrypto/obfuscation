from __future__ import print_function
from uuid import uuid4 as uuid

def straddling_sets(n):
    '''
    Creates an 'n'-straddling set system.
    '''
    symbols = [str(uuid()) for _ in xrange(2 * n - 1)]
    set1, set2 = [], []
    tmp1, tmp2 = set(), set()
    for i, symbol in enumerate(symbols, 1):
        tmp1.add(symbol)
        tmp2.add(symbol)
        if (i - 1) % 2 == 0: # i is odd
            set1.append(tmp1)
            tmp1 = set()
        else:
            set2.append(tmp2)
            tmp2 = set()
    set2.append(tmp2)
    return set1, set2

class ZimmermanObfuscator(object):
    def __init__(self, verbose=False):
        self._verbose = verbose
        _new_obf.verbose(self._verbose)
        self.logger = utils.make_logger(self._verbose)

    def obfuscate(self, circuit, secparam, directory):
        # remove old files in obfuscation directory
        if os.path.isdir(directory):
            files = os.listdir(directory)
            for file in os.listdir(directory):
                p = os.path.join(directory, file)
                os.unlink(p)

        # TODO:
        # - compute depth of 'circuit'
        # - compute deg(y)
        # - compute deg(x_i) for all i

        # - compute straddling set system

        # - construct index set
