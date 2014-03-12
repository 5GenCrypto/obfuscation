#!/usr/bin/env sage -python

from __future__ import print_function

from indobf.branchingprogram import BranchingProgram, ParseException
from indobf.obfuscator import Obfuscator

class TestParams(object):
    def __init__(self, obliviate=False, randomize=False, obfuscate=False):
        self.obliviate = obliviate
        self.randomize = randomize
        self.obfuscate = obfuscate

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
        bp = BranchingProgram(path, type='circuit', verbose=verbose)
    except ParseException as e:
        print('\x1b[33mParse Error:\x1b[0m %s' % e)
        return False
    program = bp
    if params.obliviate:
        bp.obliviate()
    if params.randomize:
        bp.randomize(secparam)
    if params.obfuscate:
        obf = Obfuscator(secparam, verbose=verbose)
        obf.obfuscate(bp)
        obf.save("%s.obf" % path)
        program = obf
    success = True
    for k, v in testcases.items():
        if program.evaluate(k) != v:
            print('\x1b[31mFail\x1b[0m (%s != %d) ' % (k, v))
            success = False
    if success:
        print('\x1b[32mPass\x1b[0m')
    return success
