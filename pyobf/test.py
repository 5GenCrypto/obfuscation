from __future__ import print_function

from pyobf.circuit import ParseException
import pyobf.utils as utils

failstr = utils.clr_error('Fail')

def test_obfuscation(path, cls, testcases, args, formula=True,
                     dual_input=False):
    success = True
    obf = cls(args.mmap, base=args.base, verbose=args.verbose,
              nthreads=args.nthreads, ncores=args.ncores)
    directory = args.save if args.save \
                else '%s.obf.%d' % (path, args.secparam)
    obf.obfuscate(path, args.secparam, directory, obliviate=args.obliviate,
                  kappa=args.kappa, formula=formula, dual_input=dual_input,
                  randomization=(not args.no_randomization))
    for k, v in testcases.items():
        if obf.evaluate(directory, k) != v:
            print('%s (%s != %d) ' % (failstr, k, v))
            success = False
    return success

def test_bp(path, cls, testcases, args):
    success = True
    try:
        if args.zimmerman:
            c = cls(path, verbose=args.verbose)
        else:
            c = cls(path, verbose=args.verbose, obliviate=args.obliviate)
    except ParseException as e:
        print('%s %s' % (utils.clr_warn('Parse Error:'), e))
        return False
    for k, v in testcases.items():
        if c.evaluate(k) != v:
            print('%s (%s != %d) ' % (failstr, k, v))
            success = False
    return success

def test_file(path, cclass, obfclass, obfuscate, args, formula=True,
              dual_input=False):
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
        success = test_obfuscation(path, obfclass, testcases, args,
                                   formula=formula, dual_input=dual_input)
    else:
        success = test_bp(path, cclass, testcases, args)
    if success:
        print(utils.clr_ok('Pass'))
    return success
