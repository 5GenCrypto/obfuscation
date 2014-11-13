from __future__ import print_function

from agis_bp import AGISBranchingProgram
from sz_bp import SZBranchingProgram
from test import test_circuit

import argparse, os, sys, time

__all__ = ['main']

errorstr = '\x1b[31mError:\x1b[0m'

def test_all(args, bpclass, obfclass, obfuscate):
    if not os.path.isdir(args.test_all):
        print("%s '%s' is not a valid directory" % (errorstr, args.test_all))
        sys.exit(1)
    ext = '.acirc' if args.zimmerman else '.circ'
    for circuit in os.listdir(args.test_all):
        path = os.path.join(args.test_all, circuit)
        if os.path.isfile(path) and path.endswith(ext):
            test_circuit(path, bpclass, obfclass, obfuscate, args)

def bp(args):
    if args.sahai_zhandry:
        cls = SZBranchingProgram
    elif args.zimmerman:
        from z_obfuscator import Circuit
        cls = Circuit
    else:
        cls = AGISBranchingProgram
    if args.test_circuit:
        test_circuit(args.test_circuit, cls, None, False, args)
    if args.test_all:
        test_all(args, cls, None, False)
    if args.load_circuit:
        bp = cls(args.load_circuit, verbose=args.verbose)
        if args.obliviate:
            bp.obliviate()
        if args.eval:
            r = bp.evaluate(args.eval)
            print('Output = %d' % r)

def obf(args):
    from agis_obfuscator import AGISObfuscator
    from sz_obfuscator import SZObfuscator
    from z_obfuscator import ZObfuscator
    if args.nslots is None:
        args.nslots = args.secparam

    if args.sahai_zhandry:
        bpclass = SZBranchingProgram
        obfclass = SZObfuscator
    elif args.zimmerman:
        bpclass = None
        obfclass = ZObfuscator
    else:
        bpclass = AGISBranchingProgram
        obfclass = AGISObfuscator

    if args.test_circuit:
        if args.attack:
            print("%s --attack flag cannot be run with --test-circuit flag" % errorstr)
            sys.exit(1)
        test_circuit(args.test_circuit, bpclass, obfclass, True, args)
    elif args.test_all:
        test_all(args, bpclass, obfclass, True)
    else:
        directory = None
        if args.load_obf:
            directory = args.load_obf
        elif args.load_circuit:
            start = time.time()
            obf = obfclass(verbose=args.verbose)
            directory = args.save if args.save \
                        else '%s.obf.%d' % (args.load_circuit, args.secparam)
            obf.obfuscate(args.load_circuit, args.secparam, directory,
                          obliviate=args.obliviate, nslots=args.nslots)
            end = time.time()
            print("Obfuscation took: %f seconds" % (end - start))
            obf.cleanup()
        else:
            print("%s One of --load-obf, --load-circuit, or --test-circuit must be used" % errorstr)
            sys.exit(1)

        if args.attack:
            assert directory
            obf = obfclass(verbose=args.verbose)
            try:
                r = obf.attack(directory, args.secparam, args.nslots)
            except AttributeError:
                print("%s --attack flag unavailable.  Make sure you compile with ATTACK = 1 in setup.py" % errorstr)
                sys.exit(1)
            print('g_1 extracted from attack: %d' % r)
        if args.eval:
            assert directory
            obf = obfclass(verbose=args.verbose)
            r = obf.evaluate(directory, args.eval)
            print('Output = %d' % r)

def main():
    parser = argparse.ArgumentParser(
        description='Cryptographic program obfuscator.')
    subparsers = parser.add_subparsers()

    parser_bp = subparsers.add_parser(
        'bp',
        help='commands for circuit -> branching program conversion')
    parser_bp.add_argument('--eval', metavar='INPUT', type=str, action='store',
                           help='evaluate branching program on INPUT')
    parser_bp.add_argument('--load-circuit', metavar='FILE', type=str,
                           action='store', help='load circuit from FILE')
    parser_bp.add_argument('--test-circuit', metavar='FILE', type=str,
                           action='store',
                           help='test BP conversion for FILE')
    parser_bp.add_argument('--test-all', metavar='DIR', nargs='?', const='circuits',
                           help='test BP conversion for all circuits in DIR (default: %(const)s)')
    parser_bp.add_argument('--secparam', metavar='N', type=int, action='store',
                           default=24, help='security parameter (default: %(default)s)')
    parser_bp.add_argument('--obliviate', action='store_true',
                           help='obliviate the branching program')
    parser_bp.add_argument('-v', '--verbose', action='store_true',
                           help='be verbose')
    parser_bp.add_argument('-s', '--sahai-zhandry', action='store_true',
                           help='use the Sahai/Zhandry construction')
    parser_bp.add_argument('-z', '--zimmerman', action='store_true',
                           help='use the Zimmerman construction')
    parser_bp.set_defaults(func=bp)

    parser_obf = subparsers.add_parser(
        'obf',
        help='commands for obfuscating a circuit/branching program')
    parser_obf.add_argument('--attack', action='store_true',
                            help='attack obfuscation')
    parser_obf.add_argument('--eval', metavar='INPUT', type=str, action='store',
                            help='evaluate obfuscation on INPUT')
    parser_obf.add_argument('--load-obf', metavar='DIR', type=str,
                            action='store',
                            help='load obfuscation from DIR')
    parser_obf.add_argument('--load-circuit', metavar='FILE', type=str,
                            action='store',
                            help='load circuit from FILE and obfuscate')
    parser_obf.add_argument('--test-circuit', metavar='FILE', type=str,
                            action='store',
                            help='test circuit from FILE')
    parser_obf.add_argument('--test-all', metavar='DIR', nargs='?', const='circuits',
                            help='test obfuscation for all circuits in DIR (default: %(const)s)')
    parser_obf.add_argument('--save', metavar='DIR', type=str, action='store',
                            help='save obfuscation to DIR')
    parser_obf.add_argument('--secparam', metavar='N', type=int, action='store',
                            default=24, help='security parameter (default: %(default)s)')
    parser_obf.add_argument('--obliviate', action='store_true',
                            help='obliviate the branching program')
    parser_obf.add_argument('--nslots', metavar='N', type=int, action='store',
                            default=None, help='number of slots to fill (default: security parameter)')
    parser_obf.add_argument('-v', '--verbose', action='store_true', 
                            help='be verbose')
    parser_obf.add_argument('-s', '--sahai-zhandry', action='store_true',
                            help='use the Sahai/Zhandry construction')
    parser_obf.add_argument('-z', '--zimmerman', action='store_true',
                            help='use the Zimmerman construction')
    parser_obf.set_defaults(func=obf)

    args = parser.parse_args()
    args.func(args)
