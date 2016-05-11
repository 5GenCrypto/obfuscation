#!/usr/bin/env python2

from __future__ import print_function
import os, random, subprocess, sys

def usage():
    print('Usage: point-json.py <bitlength> <base>')
    sys.exit(1)

def digit_to_char(digit):
    if digit < 10:
        return str(digit)
    return chr(ord('a') + digit - 10)

def str_base(number,base):
    if number < 0:
        return '-' + str_base(-number, base)
    (d, m) = divmod(number, base)
    if d > 0:
        return str_base(d, base) + digit_to_char(m)
    return digit_to_char(m)

def dary_repr(number, d, n):
    repr = list(str(str_base(number, d)))
    repr = (['0'] * (n - len(repr))) + repr
    return "".join(repr)

def random_bitstring(bitlength, base):
    num = random.randint(0, 2 ** bitlength - 1)
    repr = dary_repr(num, base, bitlength)
    bits = bin(num)[2:]
    return '0' * (bitlength - len(bits)) + bits

def run(lst):
    # print('%s' % ' '.join(lst))
    with open(os.devnull, 'w') as fnull:
        return subprocess.call(lst, stdout=fnull, stderr=fnull)

def point(bitlength, base):
    secret = random_bitstring(bitlength, base)
    lst = []
    for i in xrange(bitlength):
        lst.append('"%d"' % (i))
    str = '[' + ','.join(lst) + ']'
    fname = 'point-%d-%d.json' % (bitlength, base)
    lst = ['cryfsm', '-e', "(==) 0b%s" % secret, '-v', "\_ -> True", '-g',
           "%s" % str, '-o', fname]
    run(lst)
    with open(fname, 'r') as f:
        line = f.read()
    with open(fname, 'w') as f:
        f.write('# TEST %s 1\n' % secret)
        for _ in xrange(5):
            test = random_bitstring(bitlength, base)
            if test != secret:
                f.write('# TEST %s 0\n' % test)
        f.write(line)

def main(argv):
    if len(argv) != 3:
        usage()
    try:
        bitlength = int(argv[1])
        base = int(argv[2])
    except ValueError:
        usage()

    point(bitlength, base)

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
