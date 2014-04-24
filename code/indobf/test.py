#!/usr/bin/env sage -python

from __future__ import print_function

from branchingprogram import BranchingProgram, ParseException
from layered_branching_program import LayeredBranchingProgram
import _obfuscator as _obf

from sage.all import random_prime

def test_circuit(path, bpclass, obfclass, obfuscate, args):
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
        obf = obfclass(verbose=args.verbose, use_small_params=args.small_params,
                       use_fast_prime_gen=(not args.slow_prime_gen))
        directory = '%s.obf.%d' % (path, args.secparam)
        obf.obfuscate(bp, args.secparam, directory)
        for k, v in testcases.items():
            if obf.evaluate(directory, k) != v:
                print('\x1b[31mFail\x1b[0m (%s != %d) ' % (k, v))
                success = False
    else:
        prime = random_prime(2 ** args.secparam - 1)
        bp.randomize(prime)
        for k, v in testcases.items():
            if program.evaluate(k) != v:
                print('\x1b[31mFail\x1b[0m (%s != %d) ' % (k, v))
                success = False

    if success:
        print('\x1b[32mPass\x1b[0m')
    return success
