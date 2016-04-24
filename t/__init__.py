#!/usr/bin/env python2

from __future__ import print_function

import os, subprocess, sys

CMD = './obfuscator'
CIRCUIT_PATH = 'circuits'

yellow = '\x1b[33m'
black = '\x1b[0m'
failure_str = '\t\x1b[31mTest Failed\x1b[0m'
success_str = '\t\x1b[32mTest Succeeded\x1b[0m'

schemedict = {
    'SZ': '--sahai-zhandry',
    'Z': '--zimmerman',
}

def print_test(s):
    print('%s%s%s' % (yellow, s, black))

def run(lst):
    print('%s' % ' '.join(lst))
    return subprocess.call(lst)

def test_bp(scheme):
    print_test('Testing bp for scheme %s' % scheme)
    lst = [CMD, "bp", "--test-all", CIRCUIT_PATH]
    lst.append(schemedict[scheme])
    return run(lst)

def test_obf(scheme, mlm, secparam):
    print_test('Testing obfuscation for scheme %s' % scheme)
    lst = [CMD, "obf", "--test-all", CIRCUIT_PATH, "--secparam", str(secparam),
           "--mlm", mlm]
    lst.append(schemedict[scheme])
    return run(lst)

def test_load(scheme, mlm, secparam):
    print_test('Testing load for scheme %s' % scheme)
    if scheme == 'Z':
        circuit = 'add.acirc'
        eval = '0'
    else:
        circuit = 'and.circ'
        eval = '00'
    path = os.path.join(CIRCUIT_PATH, circuit)
    lst = [CMD, "obf", "--test", path, "--secparam", str(secparam), "--mlm", mlm]
    lst.append(schemedict[scheme])
    r = run(lst)
    if r:
        return r
    lst = [CMD, "obf", "--load-obf", path + ".obf.%d" % secparam, "--mlm", mlm, "--eval", eval]
    lst.append(schemedict[scheme])
    return run(lst)

def test(f, *args):
    if f(*args):
        print(failure_str)
    else:
        print(success_str)

def test_all():
    print("TESTING BP")
    test(test_bp, "SZ")
    test(test_bp, "Z")
    print("TESTING OBFUSCATION")
    test(test_obf, "SZ", "CLT", 12)
    test(test_obf, "SZ", "GGH", 12)
    test(test_obf, "Z", "CLT", 8)
    print("TESTING LOAD")
    test(test_load, "SZ", "CLT", 12)
    test(test_load, "SZ", "GGH", 12)
    test(test_load, "Z", "CLT", 8)

try:
    test_all()
except KeyboardInterrupt:
    pass
