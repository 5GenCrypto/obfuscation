#!/usr/bin/env sage -python

from __future__ import print_function

import argparse, sys

from indobf.test import TestParams, test_circuit

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--obliviate', action='store_true')
    args = parser.parse_args()

    paths = ['not.circ', 'and.circ', 'twoands.circ']
    params = TestParams(obliviate=args.obliviate, obfuscate=True)
    secparams = [8, 16, 24, 28, 32, 40]
    for path in paths:
        path = 'circuits/%s' % path
        print('Testing circuit "%s"' % path)
        for secparam in secparams:
            print('Security parameter = %d' % secparam)
            test_circuit(path, secparam, False, params)

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
