from __future__ import print_function

from circuit import ParseException
import zobfuscator as zobf

from sage.all import random_prime

def test_circuit(path, cclass, obfclass, obfuscate, args, zimmerman=False):
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
    success = True
    if obfuscate:
        obf = obfclass(verbose=args.verbose)
        directory = args.save if args.save \
                    else '%s.obf.%d' % (path, args.secparam)
        obf.obfuscate(path, args.secparam, directory, obliviate=args.obliviate)
        for k, v in testcases.items():
            if obf.evaluate(directory, k) != v:
                print('\x1b[31mFail\x1b[0m (%s != %d) ' % (k, v))
                success = False
    else:
        if zimmerman:
            print("\x1b[31mError:\x1b[0m Zimmerman's approach doesn't use BPs")
            return False
        # load circuit/bp
        try:
            if args.zimmerman:
                c = cclass(path, verbose=args.verbose)
            else:
                prime = random_prime(2 ** args.secparam - 1)
                c = cclass(path, verbose=args.verbose, obliviate=args.obliviate)
                c.randomize(prime)
        except ParseException as e:
            print('\x1b[33mParse Error:\x1b[0m %s' % e)
            return False
        # evaluate circuit/bp
        for k, v in testcases.items():
            # try:
                if c.evaluate(k) != v:
                    print('\x1b[31mFail\x1b[0m (%s != %d) ' % (k, v))
                    success = False
            # except Exception:
            #     print('\x1b[31mFail\x1b[0m (%s: evaluation failed)' % k)
            #     success = False
    if success:
        print('\x1b[32mPass\x1b[0m')
    return success
