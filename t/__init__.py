#!/usr/bin/env python2

from __future__ import print_function

import os, subprocess, sys

CMD = './obfuscator'
CIRCUIT_PATH = 'circuits'

yellow = '\x1b[33m'
black = '\x1b[0m'
failure_str = '\t\x1b[31mTest Failed\x1b[0m'
success_str = '\t\x1b[32mTest Succeeded\x1b[0m'

def print_test(s):
    print('%s%s%s' % (yellow, s, black))

def run(lst):
    print('%s' % ' '.join(lst))
    return subprocess.call(lst)

def test_bp():
    print_test('Testing bp')
    lst = [CMD, "bp", "--test-all", CIRCUIT_PATH]
    return run(lst)

def test_obf(mmap, secparam):
    print_test('Testing obfuscation')
    lst = [CMD, "obf", "--test-all", CIRCUIT_PATH, "--secparam", str(secparam),
           "--mmap", mmap]
    return run(lst)

def test_load(mmap, secparam):
    print_test('Testing load')
    circuit = 'and.circ'
    eval = '00'
    path = os.path.join(CIRCUIT_PATH, circuit)
    lst = [CMD, "obf", "--test", path, "--secparam", str(secparam), "--mmap", mmap]
    r = run(lst)
    if r:
        return r
    lst = [CMD, "obf", "--load-obf", path + ".obf.%d" % secparam, "--mmap", mmap, "--eval", eval]
    return run(lst)

def test(f, *args):
    if f(*args):
        print(failure_str)
    else:
        print(success_str)

def test_all():
    print("TESTING BP")
    test(test_bp)
    print("TESTING OBFUSCATION")
    test(test_obf, "CLT", 16)
    test(test_obf, "GGH", 16)
    print("TESTING LOAD")
    test(test_load, "CLT", 16)
    test(test_load, "GGH", 16)

try:
    test_all()
except KeyboardInterrupt:
    pass
