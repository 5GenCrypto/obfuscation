#!/usr/bin/env python2

from __future__ import print_function
from branchingprogram import BranchingProgram
import os, sys

testdir = 'circuits'

def test_circuit(path):
    testcases = {}
    print('Testing %s: ' % path, end='')
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
    bp = BranchingProgram(path, type='circuit')
    bp.obliviate()
    bp.randomize()
    failed = False
    for k, v in testcases.items():
        if bp.evaluate(k) != v:
            print('\x1b[31mFail\x1b[0m (%s != %d) ' % (k, v), end='')
            failed = True
    if failed:
        print()
    else:
        print('\x1b[32mPass\x1b[0m')

if __name__ == "__main__":
    if len(sys.argv) > 1:
        for path in sys.argv[1:]:
            test_circuit(path)
    else:
        for circuit in os.listdir('circuits'):
            path = os.path.join(testdir, circuit)
            test_circuit(path)
