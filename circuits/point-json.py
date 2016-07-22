#!/usr/bin/env python2

from __future__ import print_function
import math, os, random, subprocess, sys

def usage():
    print('Usage: point-json.py <length> <base>')
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

# turns a "10" in any base to a 001 000
def digit_dary_repr(num_str, d):
    L = []
    for x in num_str:
        L.append(format(int(x), 'b').zfill(int(math.ceil(math.log(d, 2)))))
    return "".join(L)

def run(lst):
    print('%s' % ' '.join(lst))
    with open(os.devnull, 'w') as fnull:
        return subprocess.call(lst, stdout=fnull, stderr=fnull)

def point(length, base):
    secret = random.randint(0, base ** length - 1)
    repr = dary_repr(secret, base, length)
    base_2_len = len(str_base(base-1, 2))
    inp = ''.join([dary_repr(int(i, base=base), 2, base_2_len) for i in repr])
    lst = []
    for i in xrange(length):
        for b in xrange(base_2_len):
            lst.append('"%d"' % i)
    str = '[' + ','.join(lst) + ']'
    fname = 'point-%d-%d.json' % (length, base)
    lst = ['cryfsm', 'util.cry', '-e', "(!=) 0b%s" % inp,
           '-v', "adjacentConstantBase `{base=%d}" % base, '-g',
           "%s" % str, '-o', fname]
    run(lst)
    with open(fname, 'r') as f:
        line = f.read()
    with open(fname, 'w') as f:
        f.write('# TEST %s 0\n' % repr)
        for _ in xrange(5):
            test = random.randint(0, base ** length - 1)
            repr = dary_repr(test, base, length)
            if test != secret:
                f.write('# TEST %s 1\n' % repr)
        f.write(line)

def main(argv):
    if len(argv) != 3:
        usage()
    try:
        length = int(argv[1])
        base = int(argv[2])
    except ValueError:
        usage()

    point(length, base)

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
