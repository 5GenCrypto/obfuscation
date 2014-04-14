#!/usr/bin/env sage -python

from __future__ import print_function

from layered_branching_program import LayeredBranchingProgram
from branchingprogram import BranchingProgram
from test import test_circuit

import argparse, os, sys, time

def bp(args):
    testdir = 'circuits'
    if args.no_layered:
        bpclass = BranchingProgram
    else:
        bpclass = LayeredBranchingProgram
    if args.test_circuit:
        test_circuit(args.test_circuit, bpclass, False, args)
    if args.test_all:
        for circuit in os.listdir('circuits'):
            path = os.path.join(testdir, circuit)
            if os.path.isfile(path) and path.endswith('.circ'):
                test_circuit(path, False, args)
    if args.load_circuit:
        bp = bpclass(args.load_circuit, verbose=args.verbose)
    if args.eval:
        r = bp.evaluate(args.eval)
        print('Output = %d' % r)

def obf(args):
    from obfuscator import Obfuscator
    if args.no_layered:
        bpclass = BranchingProgram
    else:
        bpclass = LayeredBranchingProgram
    if args.test_circuit:
        test_circuit(args.test_circuit, bpclass, True, args)
    else:
        obf = Obfuscator(verbose=args.verbose)
        if args.load_obf:
            print("Loading obfuscation from '%s'..." % args.load_obf)
            start = time.time()
            obf.load(args.load_obf)
            end = time.time()
            print("Loading took: %f seconds" % (end - start))
        elif args.load_circuit:
            print("Converting '%s' -> bp..." % args.load_circuit)
            bp = bpclass(args.load_circuit, type='circuit')
            print('Obfuscating BP of length %d...' % len(bp))
            start = time.time()
            obf.obfuscate(bp, args.secparam)
            end = time.time()
            print("Obfuscation took: %f seconds" % (end - start))
        else:
            print('One of --load-obf, --load-circuit, --test-circuit must be used!')
            sys.exit(1)

        if args.save is not None:
            print("Saving obfuscation to '%s'..." % args.save)
            obf.save(args.save)

        if args.eval is not None:
            print("Evaluating on input '%s'..." % args.eval)
            start = time.time()
            r = obf.evaluate(args.eval)
            print('Output = %d' % r)
            end = time.time()
            print("Evalution took: %f seconds" % (end - start))

def test(args):
    assert False                # XXX: not fixed up yet!
    paths = ['not.circ', 'and.circ', 'twoands.circ']
    params = TestParams(obliviate=False, obfuscate=True, fast=True)
    secparams = [8, 16, 24, 28, 32, 40]
    for path in paths:
        path = 'circuits/%s' % path
        print('Testing circuit "%s"' % path)
        for secparam in secparams:
            print('Security parameter = %d' % secparam)
            test_circuit(path, secparam, False, params)

def main():
    parser = argparse.ArgumentParser(
        description='Run indistinguishability obfuscator.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()

    parser_bp = subparsers.add_parser(
        'bp',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help='commands for circuit -> branching program conversion')
    parser_bp.add_argument('--eval', metavar='INPUT', type=str, action='store',
                           help='evaluate branching program on INPUT')
    parser_bp.add_argument('--load-circuit', metavar='FILE', type=str, action='store',
                           help='load circuit from FILE')
    parser_bp.add_argument('--test-circuit', metavar='FILE', type=str, action='store',
                           help='test FILE circuit -> bp conversion')
    parser_bp.add_argument('--test-all', action='store_true',
                           help='test circuit -> bp conversion')
    parser_bp.add_argument('--secparam', metavar='N', type=int,
                           action='store', default=8, help="security parameter")
    parser_bp.add_argument('--no-layered', action='store_true',
                           help='use standard branching programs instead of the layered variant')
    parser_bp.add_argument('-v', '--verbose', action='store_true',
                           help='be verbose')
    parser_bp.set_defaults(func=bp)

    parser_obf = subparsers.add_parser(
        'obf',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help='commands for obfuscating a circuit/branching program')
    parser_obf.add_argument('--eval', metavar='INPUT', type=str, action='store',
                            help='evaluate obfuscation on INPUT')
    parser_obf.add_argument('--load-obf', metavar='DIR', type=str, action='store',
                            help='load obfuscation from DIR')
    parser_obf.add_argument('--load-circuit', metavar='FILE', type=str, action='store',
                            help='load circuit from FILE and obfuscate')
    parser_obf.add_argument('--test-circuit', metavar='FILE', type=str, action='store',
                            help='test FILE -> obfuscation')
    parser_obf.add_argument('--save', metavar='DIR', type=str, action='store',
                            help='save obfuscation to DIR')
    parser_obf.add_argument('--secparam', metavar='N', type=int,
                             action='store', default=8, help="security parameter")
    parser_obf.add_argument('--fast', action='store_true',
                             help='use smaller parameters so things run faster')
    parser_obf.add_argument('-v', '--verbose', action='store_true',
                            help='be verbose')
    parser_obf.set_defaults(func=obf)

    parser_test = subparsers.add_parser(
        'test',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help='run test suite')
    parser_test.set_defaults(func=test)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        pass
