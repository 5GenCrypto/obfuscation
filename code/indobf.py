#!/usr/bin/env sage -python

from __future__ import print_function

from indobf.branchingprogram import BranchingProgram
from indobf.gradedencoding import GradedEncoding
from indobf.obfuscator import Obfuscator
from indobf.test import TestParams, test_circuit

from sage.all import *

import argparse, os, sys, time

def bp(args):
    testdir = 'circuits'
    params = TestParams(obliviate=True)
    if args.test_circuit:
        test_circuit(args.test_circuit, args.secparam, args.verbose, params)
    if args.test_all:
        for circuit in os.listdir('circuits'):
            path = os.path.join(testdir, circuit)
            test_circuit(path, args.secparam, args.verbose, params)
    if args.load_circuit:
        bp = BranchingProgram(args.load_circuit, type='circuit')
    if args.eval:
        r = bp.evaluate(args.eval)
        print('Output = %d' % r)

def obf(args):
    if args.test_circuit:
        params = TestParams(obliviate=False, obfuscate=True,
                            disable_mbundling=args.disable_mbundling,
                            disable_bookends=args.disable_bookends)
        test_circuit(args.test_circuit, args.secparam, args.verbose, params)
    else:
        obf = Obfuscator(verbose=args.verbose)
        if args.load_obf is not None:
            print("Loading obfuscation from '%s'..." % args.load_obf)
            start = time.time()
            obf.load(args.load_obf)
            end = time.time()
            print("Loading took: %f seconds" % (end - start))
        elif args.load_circuit is not None:
            print("Converting '%s' -> bp..." % args.load_circuit)
            bp = BranchingProgram(args.load_circuit, type='circuit')
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

def main(argv):
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
    parser_obf.add_argument('--disable-mbundling', action='store_true',
                             help='disable multiplicative bundling')
    parser_obf.add_argument('--disable-bookends', action='store_true',
                             help='disable booken vectors')
    parser_obf.add_argument('-v', '--verbose', action='store_true',
                            help='be verbose')
    parser_obf.set_defaults(func=obf)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
