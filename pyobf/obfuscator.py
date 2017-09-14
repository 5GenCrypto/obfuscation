import pyobf._obfuscator as _obf
from pyobf.sz_bp import SZBranchingProgram
import pyobf.utils as utils
import os, re, time

MMAP_CLT = 0x00
MMAP_GGHLITE = 0x01
MMAP_DUMMY = 0x02

OBFUSCATOR_FLAG_NONE = 0x00
OBFUSCATOR_FLAG_NO_RANDOMIZATION = 0x01
OBFUSCATOR_FLAG_VERBOSE = 0x04

ENCODE_LAYER_RANDOMIZATION_TYPE_NONE = 0x00
ENCODE_LAYER_RANDOMIZATION_TYPE_FIRST = 0x01
ENCODE_LAYER_RANDOMIZATION_TYPE_MIDDLE = 0x02
ENCODE_LAYER_RANDOMIZATION_TYPE_LAST = 0x04

def get_mmap_flag(mmap):
    if mmap == 'CLT':
        return MMAP_CLT
    elif mmap == 'GGH':
        return MMAP_GGHLITE
    else:
        return MMAP_DUMMY

class Obfuscator(object):
    def __init__(self, mmap, base=None, verbose=False, nthreads=None,
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

    '''
    Get size of obfuscation (in bytes)
    '''
    def obfsize(self, directory):
        size = 0
        for f in os.listdir(directory):
            size += os.path.getsize(os.path.join(directory, f))
        return size

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
        self.logger('Obfuscation took: %f s' % (end - start))
        self.logger('Obfuscation size: %0.2f KB' % (self.obfsize(directory) / 1024.0))
        if self._verbose:
            _obf.max_mem_usage()

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
        files = os.listdir(directory)
        if self._base:
            base = self._base
        else:
            # Compute base by counting the number of files that correspond to
            # the first MBP layer.
            base = len(filter(lambda s: re.match('0.\d+', s), files))
        # Input length is equal to the number of `[num].input` files in
        # `directory`.
        inplen = len(filter(lambda file: re.match('\d+.input', file), files))
        if len(inp) != inplen:
            print('{} Invalid input length ({} != {})'.format(
                utils.clr_error('Error:'), len(inp), inplen))
            return None
        try:
            inp = [int(i, base) for i in inp]
        except ValueError:
            print('{} Invalid input for base {}'.format(
                utils.clr_error('Error:'), base))
            return None
        flags = OBFUSCATOR_FLAG_NONE
        if self._verbose:
            flags |= OBFUSCATOR_FLAG_VERBOSE
        return self._evaluate(directory, inp, _obf.evaluate, _obf, flags)
