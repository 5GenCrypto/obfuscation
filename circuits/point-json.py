#!/usr/bin/env python2

from __future__ import print_function
import  argparse, os, random, sys
from util import *

def point(base, length, cryfsm='cryfsm'):
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
    lst = [cryfsm, 'util.cry', '-e', "(!=) 0b%s" % inp,
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
    parser = argparse.ArgumentParser()
    parser.add_argument('--cryfsm', metavar='PATH', action='store', type=str,
                        default='cryfsm',
                        help='set PATH as cryfsm executable')
    parser.add_argument('base', action='store', type=int)
    parser.add_argument('length', action='store', type=int)
    args = parser.parse_args()
    point(args.base, args.length, args.cryfsm)

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
