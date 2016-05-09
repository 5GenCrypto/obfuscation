#!/usr/bin/env python2

from __future__ import print_function
import os, random, subprocess, sys

def usage():
    print('Usage: point-json.py <bitlength>')
    sys.exit(1)

def random_bitstring(bitlength):
    r = random.randint(0, 2 ** bitlength - 1)
    bits = bin(r)[2:]
    return '0' * (bitlength - len(bits)) + bits

def run(lst):
    # print('%s' % ' '.join(lst))
    with open(os.devnull, 'w') as fnull:
        return subprocess.call(lst, stdout=fnull, stderr=fnull)

def binary_point(bitlength):
    secret = random_bitstring(bitlength)
    lst = ['cryfsm', '-e', "(==) 0b%s" % secret, '-v', "\_ -> True", '-g',
           '#', '-o', 'point-%d.json' % bitlength]
    run(lst)
    with open('point-%d.json' % bitlength, 'r') as f:
        line = f.read()
    with open('point-%d.json' % bitlength, 'w') as f:
        f.write('# TEST %s 1\n' % secret)
        for _ in xrange(5):
            test = random_bitstring(bitlength)
            if test != secret:
                f.write('# TEST %s 0\n' % test)
        f.write(line)

def main(argv):
    if len(argv) != 2:
        usage()
    try:
        bitlength = int(argv[1])
    except ValueError:
        usage()

    binary_point(bitlength)

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
