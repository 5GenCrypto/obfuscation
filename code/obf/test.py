from __future__ import print_function

from circuit import ParseException
from sage.all import random_prime

__all__ = ['test_circuit']

failstr = '\x1b[31mFail\x1b[0m'
parseerrstr = '\x1b[33mParse Error:\x1b[0m'

def test_obfuscation(path, cls, testcases, args):
    success = True
    obf = cls(verbose=args.verbose)
    directory = args.save if args.save \
                else '%s.obf.%d' % (path, args.secparam)
    obf.obfuscate(path, args.secparam, directory, obliviate=args.obliviate,
                  nslots=args.nslots, kappa=args.kappa)
    for k, v in testcases.items():
        r = obf.evaluate(directory, k)
        if r != v:
            print('%s (%s != %d) ' % (failstr, k, v))
            success = False
    return success

def test_bp(path, cls, testcases, args):
    success = True
    try:
        if args.zimmerman:
            c = cls(path, verbose=args.verbose)
        else:
            prime = random_prime(2 ** args.secparam - 1)
            c = cls(path, verbose=args.verbose, obliviate=args.obliviate)
            c.randomize(prime)
    except ParseException as e:
        print('%s %s' % (parseerrstr, e))
        return False
    for k, v in testcases.items():
        if c.evaluate(k) != v:
            print('%s (%s != %d) ' % (failstr, k, v))
            success = False
    return success

def test_circuit(path, cclass, obfclass, obfuscate, args):
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
    if obfuscate:
        success = test_obfuscation(path, obfclass, testcases, args)
    else:
        success = test_bp(path, cclass, testcases, args)
    if success:
        print('\x1b[32mPass\x1b[0m')
    return success
