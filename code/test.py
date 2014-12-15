#!/usr/bin/env python2

from __future__ import print_function

import os, subprocess

CMD = './obfuscator'
CIRCUIT_PATH = 'circuits'

yellow = '\x1b[33m'
black = '\x1b[0m'
failure_str = '\t\x1b[31mTest Failed\x1b[0m'
success_str = '\t\x1b[32mTest Succeeded\x1b[0m'

schemedict = {
    'AGIS': '',
    'SZ': '--sahai-zhandry',
    'Z': '--zimmerman',
}

def print_test(s):
    print('%s%s%s' % (yellow, s, black))

def run(lst):
    print('%s' % ' '.join(lst))
    return subprocess.call(lst)

def test_circuits(cmd, scheme):
    print_test('Testing all circuits for scheme %s' % scheme)
    lst = [cmd, "obf", "--test-all", CIRCUIT_PATH, "--secparam", "8"]
    if scheme != 'AGIS':
        lst.append(schemedict[scheme])
    return run(lst)

def test_load(cmd, scheme):
    print_test('Testing load for scheme %s' % scheme)
    if scheme == 'Z':
        circuit = 'add.acirc'
        eval = '0'
    else:
        circuit = 'and.circ'
        eval = '00'
    path = os.path.join(CIRCUIT_PATH, circuit)
    lst = [cmd, "obf", "--test-circuit", path, "--secparam", "8"]
    if scheme != 'AGIS':
        lst.append(schemedict[scheme])
    r = run(lst)
    if r:
        return r
    lst = [cmd, "obf", "--load-obf", path + ".obf.8", "--eval", eval]
    if scheme != 'AGIS':
        lst.append(schemedict[scheme])
    return run(lst)

def test(f, cmd, *args):
    if f(cmd, *args):
        print(failure_str)
    else:
        print(success_str)

def test_all(cmd):
    test(test_load, cmd, "AGIS")
    test(test_load, cmd, "SZ")
    test(test_load, cmd, "Z")
    test(test_circuits, cmd, "AGIS")
    test(test_circuits, cmd, "SZ")
    test(test_circuits, cmd, "Z")

if __name__ == '__main__':
    try:
        test_all(CMD)
    except KeyboardInterrupt:
        pass
