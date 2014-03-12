#!/usr/bin/env sage -python

from __future__ import print_function

from indobf.test import TestParams, test_circuit

def main():
    paths = ['circuits/not.circuit', 'circuits/and.circuit']
    params = TestParams(obfuscate=True)
    secparams = [8, 12, 16, 20, 24, 28]
    for path in paths:
        print('Testing circuit "%s"' % path)
        for secparam in secparams:
            print('Security parameter = %d' % secparam)
            test_circuit(path, secparam, False, params)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
