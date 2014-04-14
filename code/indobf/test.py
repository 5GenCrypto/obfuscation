#!/usr/bin/env sage -python

from __future__ import print_function

from branchingprogram import BranchingProgram, ParseException
from layered_branching_program import LayeredBranchingProgram
import _obfuscator as _obf

def test_circuit(path, bpclass, obfuscate, args):
    testcases = {}
    print('Testing %s: ' % path, end='')
    if args.verbose:
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
        bp = bpclass(path, verbose=args.verbose)
    except ParseException as e:
        print('\x1b[33mParse Error:\x1b[0m %s' % e)
        return False
    program = bp
    # if args.obliviate:
    #     bp.obliviate()
    success = True
    if obfuscate:
        from obfuscator import Obfuscator
        kwargs = {
            'verbose': args.verbose,
            'fast': args.fast,
        }
        obf = Obfuscator(**kwargs)
        directory = '%s.obf' % path
        obf.obfuscate(bp, args.secparam, directory)
        for k, v in testcases.items():
            if obf.evaluate(directory, k) != v:
                print('\x1b[31mFail\x1b[0m (%s != %d) ' % (k, v))
                success = False
    else:
        # prime = _obf.genprime(args.secparam)
        # bp.randomize(prime, alphas=None)
        for k, v in testcases.items():
            if program.evaluate(k) != v:
                print('\x1b[31mFail\x1b[0m (%s != %d) ' % (k, v))
                success = False

    if success:
        print('\x1b[32mPass\x1b[0m')
    return success
