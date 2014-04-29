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
        bp = bpclass(path, verbose=args.verbose, obliviate=args.obliviate)
    except ParseException as e:
        print('\x1b[33mParse Error:\x1b[0m %s' % e)
        return False
    program = bp
    success = True
    if obfuscate:
        def create_obf():
            return obfclass(verbose=args.verbose,
                            use_fast_prime_gen=(not args.slow_prime_gen))
        obf = create_obf()
        directory = args.save if args.save \
                    else '%s.obf.%d' % (path, args.secparam)
        obf.obfuscate(bp, args.secparam, directory)
        obf.encode_benchmark()
        obf.cleanup()
        for k, v in testcases.items():
            if obf.evaluate(directory, k) != v:
                print('\x1b[31mFail\x1b[0m (%s != %d) ' % (k, v))
                success = False
    else:
        prime = random_prime(2 ** args.secparam - 1)
        bp.randomize(prime)
        for k, v in testcases.items():
            try:
                if program.evaluate(k) != v:
                    print('\x1b[31mFail\x1b[0m (%s != %d) ' % (k, v))
                    success = False
            except Exception:
                print('\x1b[31mFail\x1b[0m (%s: evaluation failed)' % k)
                success = False

    if success:
        print('\x1b[32mPass\x1b[0m')
    return success
