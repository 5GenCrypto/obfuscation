from __future__ import print_function

from pyobf.circuit import ParseException
from pyobf.test import test_file
from pyobf.sz_bp import SZBranchingProgram
from pyobf.obfuscator import Obfuscator

import argparse, os, sys, time
import pyobf.utils as utils

__all__ = ['main']

errorstr = utils.clr_error('Error:')

def is_formula(fname, args):
    ext = os.path.splitext(fname)[1]
    if ext in ('.circ'):
        return True
    elif ext == '.json':
        return False
    else:
        print("%s unknown extension '%s'" % (errorstr, ext))
        sys.exit(1)

def test_all(args, obfuscate):
    success = True
    if not os.path.isdir(args.test_all):
        print("%s '%s' is not a valid directory" % (errorstr, args.test_all))
        sys.exit(1)
    ext = '.circ'
    for circuit in os.listdir(args.test_all):
        path = os.path.join(args.test_all, circuit)
        if os.path.isfile(path) and path.endswith(ext):
            success &= test_file(path, obfuscate, args)
        if os.path.isfile(path) and path.endswith('.json'):
            success &= test_file(path, obfuscate, args, formula=False)
    return success

def bp(args):
    success = True
    try:
        if args.test:
            success = test_file(args.test, False, args)
        elif args.test_all:
            success = test_all(args, False)
        elif args.load:
            formula = is_formula(args.load, args)
            bp = SZBranchingProgram(args.load, base=args.base,
                                    verbose=args.verbose, formula=formula)
            if args.print:
                print(bp)
            if args.eval:
                r = bp.evaluate(args.eval)
                print('Output = %d' % r)
    except ParseException as e:
        print('%s %s' % (errorstr, e))
        sys.exit(1)
    return success

def obf(args):
    if args.mmap not in ('CLT', 'GGH', 'DUMMY'):
        print('--mmap must be either CLT, GGH, or DUMMY')
        sys.exit(1)
    success = True
    try:
        if args.test:
            formula = is_formula(args.test, args)
            success = test_file(args.test, True, args, formula=formula)
        elif args.test_all:
            success = test_all(args, True)
        else:
            directory = None
            if args.load_obf:
                directory = args.load_obf
            elif args.load:
                formula = is_formula(args.load, args)
                obf = Obfuscator(args.mmap, base=args.base,
                                 verbose=args.verbose, nthreads=args.nthreads,
                                 ncores=args.ncores)
                directory = args.save if args.save \
                            else '%s.obf.%d' % (args.load, args.secparam)
                obf.obfuscate(args.load, args.secparam, directory,
                              kappa=args.kappa, formula=formula,
                              randomization=(not args.no_randomization),
                              seed=args.seed)
            else:
                print('%s One of --load-obf, --load, or '
                      '--test must be used' % errorstr)
                sys.exit(1)

            if args.eval:
                assert directory
                obf = Obfuscator(args.mmap, base=args.base,
                                 verbose=args.verbose, nthreads=args.nthreads,
                                 ncores=args.ncores)
                r = obf.evaluate(directory, args.eval)
                if r is not None:
                    print('Output = %d' % r)
    except ParseException as e:
        print('%s %s' % (errorstr, e))
        sys.exit(1)
    return success

def main():
    parser = argparse.ArgumentParser(
        description='Cryptographic program obfuscator.')
    subparsers = parser.add_subparsers()

    try:
        ncores = os.sysconf('SC_NPROCESSORS_ONLN')
    except ValueError:
        print(utils.clr_warn('Warning: Unable to count number of cores, defaulting to 1'))
        ncores = 1
    nthreads = ncores
    secparam = 24

    parser_bp = subparsers.add_parser(
        'bp',
        help='commands for circuit -> branching program conversion')
    parser_bp.add_argument('--eval',
                           metavar='INPUT', action='store', type=str,
                           help='evaluate branching program on INPUT')
    parser_bp.add_argument('--load',
                           metavar='FILE', action='store', type=str,
                           help='load circuit or branching program from FILE')
    parser_bp.add_argument('--test',
                           metavar='FILE', action='store', type=str,
                           help='test branching program conversion for FILE')
    parser_bp.add_argument('--test-all',
                           metavar='DIR', nargs='?', const='circuits/',
                           help='test branching program conversion for all circuits in DIR (default: %(const)s)')
    parser_bp.add_argument('--base',
                           metavar='B', action='store', type=int, default=None,
                           help='base of matrix branching program (default: guess)')
    parser_bp.add_argument('--print',
                           action='store_true',
                           help='print branching program to stdout')
    parser_bp.add_argument('-v', '--verbose',
                           action='store_true',
                           help='be verbose')
    parser_bp.set_defaults(func=bp)

    parser_obf = subparsers.add_parser(
        'obf',
        help='commands for obfuscating a circuit/branching program')
    parser_obf.add_argument('--eval',
                            metavar='INPUT', action='store', type=str,
                            help='evaluate obfuscation on INPUT')
    parser_obf.add_argument('--kappa',
                            metavar='N', action='store', type=int, default=None,
                            help='set kappa to N (for debugging)')
    parser_obf.add_argument('--load-obf',
                            metavar='DIR', action='store', type=str,
                            help='load obfuscation from DIR')
    parser_obf.add_argument('--load',
                            metavar='FILE', action='store', type=str,
                            help='load circuit or branching program from FILE')
    parser_obf.add_argument('--test',
                            metavar='FILE', action='store', type=str,
                            help='test circuit or branching program from FILE')
    parser_obf.add_argument('--test-all',
                            metavar='DIR', nargs='?', const='circuits/',
                            help='test obfuscation for all circuits in DIR (default: %(const)s)')
    parser_obf.add_argument('--mmap', metavar='M', type=str, default='CLT',
                            action='store',
                            help='use multilinear map M [CLT, GGH, DUMMY] (default: %(default)s)')
    parser_obf.add_argument('--save',
                            metavar='DIR', action='store', type=str,
                            help='save obfuscation to DIR')
    parser_obf.add_argument('--secparam',
                            metavar='N', action='store', type=int,
                            default=secparam, help='security parameter (default: %(default)s)')
    parser_obf.add_argument('--nthreads',
                            metavar='N', action='store', type=int, default=nthreads,
                            help='number of threads to use in threadpool (default: %(default)s)')
    parser_obf.add_argument('--ncores',
                            metavar='N', action='store', type=int, default=ncores,
                            help='number of cores to use for OpenMP (default: %(default)s)')
    parser_obf.add_argument('--base',
                            metavar='B', action='store', type=int, default=None,
                            help='base of matrix branching program (default: guess)')
    parser_obf.add_argument('--seed',
                            metavar='FILE', action='store', type=str,
                            help='load seed from FILE')
    parser_obf.add_argument('--no-randomization', action='store_true',
                            help='turn of branching program randomization')
    parser_obf.add_argument('-v', '--verbose',
                            action='store_true',
                            help='be verbose')
    parser_obf.set_defaults(func=obf)

    args = parser.parse_args()
    return args.func(args)
