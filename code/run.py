#!/usr/bin/env sage -python

from __future__ import print_function

from branchingprogram import BranchingProgram
from obfuscator import Obfuscator

from sage.all import *

import argparse, sys, time


def main(argv):
    parser = argparse.ArgumentParser(
        description='Run obfuscator.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('circuit', type=str)
    parser.add_argument('input', type=str)
    parser.add_argument('--ncpus', metavar='N', type=int, action='store',
                        default=sage.parallel.ncpus.ncpus(),
                        help='set number of CPUs to use (only used when --parallel flag set)')
    parser.add_argument('-p', '--parallel', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('--secparam', metavar='N', type=int, action='store',
                        default=8, help="set security parameter")
    args = parser.parse_args()

    start = time.time()
    print('Converting circuit -> bp...')
    bp = BranchingProgram(args.circuit, type='circuit')
    # bp.obliviate()
    bp.randomize()
    print('Obfuscating BP of length %d...' % len(bp))
    obf = Obfuscator(args.secparam, bp, verbose=args.verbose,
                     parallel=args.parallel, ncpus=args.ncpus)
    obf.obfuscate()
    print('Evaluating on input %s...' % args.input)
    r = obf.evaluate(args.input)
    print('Output = %d' % r)
    end = time.time()
    print('Total time: %f seconds' % (end - start))


if __name__ == "__main__":
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
