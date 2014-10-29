from __future__ import print_function

import _obfuscator as _obf
import utils
from bp_sz import SZBranchingProgram

import os, time

def to_long(lst):
    return [long(i) for i in lst]

def pad(array, length, bplength):
    if len(array) < length:
        zeros = [to_long([0] * bplength)]
        return array + (zeros * (length - len(array)))
    else:
        return array

class SZObfuscator(object):

    def __init__(self, verbose=False):
        self.obfuscation = None
        self._verbose = verbose
        _obf.verbose(self._verbose)
        self.logger = utils.make_logger(self._verbose)

    def _gen_mlm_params(self, secparam, kappa, nzs, directory):
        self.logger('Generating MLM parameters...')
        start = time.time()
        if not os.path.exists(directory):
            os.mkdir(directory)
        self._state, primes = _obf.setup(secparam, kappa, 0, nzs, directory)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return primes

    def _construct_bps(self, nslots, circuit, primes, obliviate):
        self.logger('Constructing %d BP...' % nslots)
        start = time.time()
        bps = []
        for _, prime in zip(xrange(nslots), primes):
            bp = SZBranchingProgram(circuit, verbose=self._verbose,
                                    obliviate=obliviate)
            bp.set_straddling_sets()
            bps.append(bp)
        end = time.time()
        self.logger('Took: %f' % (end - start))
        return bps

    def _randomize(self, secparam, bps, primes):
        self.logger('Randomizing BPs...')
        start = time.time()
        for bp, prime in zip(bps, primes):
            # bp.randomize(prime)
            bp.set_straddling_sets()
        end = time.time()
        self.logger('Took: %f' % (end - start))

    def _obfuscate(self, bps, length):
        print(bps[0])
        print(len(bps[0]))
        for i in xrange(len(bps[0])):
            self.logger('Obfuscating layer...')
            start = time.time()
            bplength = len(bps[0][i].zero.transpose().list())
            zeros = pad([to_long(bp[i].zero.transpose().list()) for bp in bps],
                        length, bplength)
            ones = pad([to_long(bp[i].one.transpose().list()) for bp in bps],
                       length, bplength)
            nrows = bps[0][i].zero.nrows()
            ncols = bps[0][i].zero.ncols()
            _obf.encode_layers(self._state, i, nrows, ncols, bps[0][i].inp,
                               zeros, ones, bps[0][i].zeroset, bps[0][i].oneset)
            end = time.time()
            self.logger('Took: %f' % (end - start))
    
    def obfuscate(self, circuit, secparam, directory, obliviate=False,
                  nslots=None):
        # remove old files in obfuscation directory
        if os.path.isdir(directory):
            files = os.listdir(directory)
            for file in os.listdir(directory):
                p = os.path.join(directory, file)
                os.unlink(p)

        if nslots is None:
            nslots = secparam

        # create a dummy branching program to determine parameters
        bp = SZBranchingProgram(circuit, verbose=self._verbose,
                                obliviate=obliviate)
        kappa = len(bp)
        nzs = bp.set_straddling_sets()

        start = time.time()
        primes = self._gen_mlm_params(secparam, kappa, nzs, directory)
        bps = self._construct_bps(nslots, circuit, primes, obliviate)
        self._randomize(secparam, bps, primes)
        self._obfuscate(bps, len(primes))
        end = time.time()
        self.logger('Obfuscation took: %f' % (end - start))
        if self._verbose:
            _obf.max_mem_usage()

    def evaluate(self, directory, inp):
        self.logger('Evaluating %s...' % inp)
        start = time.time()
        files = os.listdir(directory)
        inputs = sorted(filter(lambda s: 'input' in s, files))
        result = _obf.sz_evaluate(directory, inp, len(inputs))
        end = time.time()
        self.logger('Took: %f' % (end - start))
        if self._verbose:
            _obf.max_mem_usage()
        return result

    def cleanup(self):
        _obf.cleanup(self._state)
