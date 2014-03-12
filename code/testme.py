#!/usr/bin/env sage -python

from __future__ import print_function

import argparse, sys

from indobf.test import TestParams, test_circuit

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--obliviate', action='store_true')
    parser.add_argument('--randomize', action='store_true')
    args = parser.parse_args()

    paths = ['circuits/not.circuit', 'circuits/and.circuit']
    params = TestParams(obliviate=args.obliviate, randomize=args.randomize,
                        obfuscate=True)
    secparams = [8, 12, 16, 20, 24, 28]
    for path in paths:
        print('Testing circuit "%s"' % path)
        for secparam in secparams:
            print('Security parameter = %d' % secparam)
            test_circuit(path, secparam, False, params)

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
