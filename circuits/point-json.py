#!/usr/bin/env python2

from __future__ import print_function
import  os, random, sys
from util import *

def usage(ret):
    print('Usage: point-json.py <base> <length>')
    sys.exit(ret)

def point(base, length):
    secret = random.randint(0, base ** length - 1)
    repr = dary_repr(secret, base, length)
    base_2_len = len(str_base(base-1, 2))
    inp = ''.join([dary_repr(int(i, base=base), 2, base_2_len) for i in repr])
    lst = []
    for i in xrange(length):
        for b in xrange(base_2_len):
            lst.append('"%d"' % i)
    str = '[' + ','.join(lst) + ']'
    tmp = 'point-%d-%d-tmp.json' % (base, length)
    fname = 'point-%d-%d.json' % (base, length)
    lst = ['cryfsm', 'util.cry', '-e', "(!=) 0b%s" % inp,
           '-v', "adjacentConstantBase `{base=%d}" % base, '-g',
           "%s" % str, '-o', tmp]
    run(lst)
    lst = ['fsmevade', '-i', tmp, '-o', fname, 'True']
    run(lst)
    os.remove(tmp)
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
        usage(1)
    try:
        base = int(argv[1])
        length = int(argv[2])
    except ValueError:
        usage(1)

    point(base, length)

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
