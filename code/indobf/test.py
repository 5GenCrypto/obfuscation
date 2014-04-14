#!/usr/bin/env sage -python

from __future__ import print_function

from branchingprogram import BranchingProgram, ParseException
from layered_branching_program import LayeredBranchingProgram
import _obfuscator as _obf

class TestParams(object):
    def __init__(self, obliviate=False, obfuscate=False, fast=False):
        self.obliviate = obliviate
        self.obfuscate = obfuscate
        self.fast = fast

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
        bp = LayeredBranchingProgram(path, verbose=verbose)
    except ParseException as e:
        print('\x1b[33mParse Error:\x1b[0m %s' % e)
        return False
    program = bp
    if params.obliviate:
        bp.obliviate()
    success = True
    if params.obfuscate:
        from obfuscator import Obfuscator
        kwargs = {
            'verbose': verbose,
            'fast': params.fast,
        }
        obf = Obfuscator(**kwargs)
        directory = '%s.obf' % path
        obf.obfuscate(bp, secparam, directory)
        for k, v in testcases.items():
            if obf.evaluate(directory, k) != v:
                print('\x1b[31mFail\x1b[0m (%s != %d) ' % (k, v))
                success = False
    else:
        # prime = _obf.genprime(secparam)
        # bp.randomize(prime, alphas=None)
        for k, v in testcases.items():
            if program.evaluate(k) != v:
                print('\x1b[31mFail\x1b[0m (%s != %d) ' % (k, v))
                success = False

    if success:
        print('\x1b[32mPass\x1b[0m')
    return success
