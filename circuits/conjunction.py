#!/usr/bin/env python2

from __future__ import print_function
import random, re, sys
from util import *

def usage(ret):
    print("""
Usage: conjunction.py <bitstring>

Produces the matrix branching program for conjunction on a given bitstring,
where 'bitstring' must contain only characters 0, 1, and ?.
Example: conjunction.py 01??1?11101?01""")
    sys.exit(ret)

def main(argv):
    if len(argv) != 2:
        usage(1)
    bitstring = argv[1]
    r = re.search("[^01?]", bitstring)
    if r:
        print("Error: invalid character '%c' in bitstring" % r.group())
        usage(1)
    length = len(bitstring)
    fname = '%s.json' % re.escape(bitstring).replace('\\', '')
    lst = ['cryfsm', 'conjunction.cry', '-g', '#', '-e', 'main \"%s\"' % bitstring,
           '-o', fname]
    run(lst)
    with open(fname, 'r') as f:
        line = f.read()
    pattern = bitstring.replace('?', '.?')
    with open(fname, 'w') as f:
        f.write('# TEST %s 1\n' % bitstring.replace('?', '1'))
        for _ in xrange(5):
            test = random.randint(0, 2 ** length - 1)
            string = dary_repr(test, 2, length)
            match = re.match(pattern, string)
            if match is None or len(match.group()) != length:
                result = 0
            else:
                result = 1
            f.write('# TEST %s %d\n' % (string, result))
        f.write(line)

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
