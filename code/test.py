#!/usr/bin/env python2

from __future__ import print_function
import branchingprogram
import os

testdir = 'circuits'

if __name__ == "__main__":
    for circuit in os.listdir('circuits'):
        path = os.path.join(testdir, circuit)
        testcases = {}
        print('Testing %s: ' % circuit, end='')
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
            continue
        bp = branchingprogram.BranchingProgram(path, type='circuit')
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
