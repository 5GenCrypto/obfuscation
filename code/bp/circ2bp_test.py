#!/usr/bin/env python2

from __future__ import print_function
import circ2bp
import os

testdir = 'circuits'

if __name__ == "__main__":
    for circuit in os.listdir('circuits'):
        path = os.path.join(testdir, circuit)
        nins = circ2bp.n_inputs(path)
        depth = None
        testcases = {}
        print('Testing %s: ' % circuit, end='')
        with open(path) as f:
            for line in f:
                if line.startswith('#'):
                    if 'DEPTH' in line:
                        _, _, depth = line.split()
                        depth = int(depth)
                    elif 'TEST' in line:
                        _, _, inp, outp = line.split()
                        testcases[inp] = int(outp)
                else:
                    continue
        if len(testcases) == 0 or nins is None or depth is None:
            print()
            continue
        bp = circ2bp.circuit_to_bp(path)
        bp = circ2bp.obliviate(bp, nins, depth)
        bp = circ2bp.randomize(bp)
        failed = False
        for k, v in testcases.items():
            if circ2bp.eval_bp(bp, k, circ2bp.MSZp) != v:
                print('\x1b[31mFail\x1b[0m (%s != %d) ' % (k, v), end='')
                failed = True
        if failed:
            print()
        else:
            print('\x1b[32mPass\x1b[0m')
