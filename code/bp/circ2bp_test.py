#!/usr/bin/env python3

import circ2bp
import os

testdir = 'circuits'

if __name__ == "__main__":
    for circuit in os.listdir('circuits'):
        nins = None
        testcases = {}
        path = os.path.join(testdir, circuit)
        print('Testing %s: ' % circuit, end='')
        with open(path) as f:
            for line in f:
                if line.startswith('#'):
                    if 'NINS' in line:
                        _, _, nins = line.split()
                        nins = int(nins)
                    elif 'TEST' in line:
                        _, _, inp, outp = line.split()
                        testcases[inp] = int(outp)
                else:
                    continue
        if len(testcases) == 0 or nins is None:
            print()
            continue
        bp = circ2bp.circuit_to_bp(path)
        obp = circ2bp.obliviate(bp, nins)
        failed = False
        for k, v in testcases.items():
            if circ2bp.eval_bp(obp, k) != v:
                print('\x1b[31m✘\x1b[0m (%s != %d) ' % (k, v), end='')
                failed = True
        if failed:
            print()
        else:
            print('\x1b[32m✔\x1b[0m')
