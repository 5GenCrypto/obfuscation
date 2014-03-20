#!/usr/bin/env sage -python

from __future__ import print_function

from indobf.branchingprogram import BranchingProgram, ParseException
from sage.rings.arith import random_prime

class TestParams(object):
    def __init__(self, obliviate=False, obfuscate=False,
                 disable_mbundling=False, disable_bookends=False):
        self.obliviate = obliviate
        self.obfuscate = obfuscate
        self.disable_mbundling = disable_mbundling
        self.disable_bookends = disable_bookends

def test_circuit(path, secparam, verbose, params):
    testcases = {}
    print('Testing %s: ' % path, end='')
    if verbose:
        print()
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                if 'TEST' in line:
                    _, _, inp, outp = line.split()
                    testcases[inp] = int(outp)
            else:
                continue
    if len(testcases) == 0:
        print('no test cases')
        return
    try:
        bp = BranchingProgram(path, type='circuit', group='S6',
                              verbose=verbose)
    except ParseException as e:
        print('\x1b[33mParse Error:\x1b[0m %s' % e)
        return False
    program = bp
    if params.obliviate:
        bp.obliviate()
    if params.obfuscate:
        from indobf.obfuscator import Obfuscator
        kwargs = {
            'verbose': verbose,
            'disable_mbundling': params.disable_mbundling,
            'disable_bookends': params.disable_bookends
        }
        obf = Obfuscator(**kwargs)
        obf.obfuscate(bp, secparam)
        obf.save('%s.obf' % path)
        program = Obfuscator(**kwargs)

        program.load('%s.obf' % path)
    else:
        prime = long(random_prime((1 << secparam) - 1, lbound=(1 << secparam - 1)))
        bp.randomize(prime, alphas=None)
    success = True
    for k, v in testcases.items():
        if program.evaluate(k) != v:
            print('\x1b[31mFail\x1b[0m (%s != %d) ' % (k, v))
            success = False
    if success:
        print('\x1b[32mPass\x1b[0m')
    return success
