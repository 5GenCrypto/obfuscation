#!/usr/bin/env python2

from __future__ import print_function
import argparse, random, re, sys
from util import *

def main(argv):
    parser = argparse.ArgumentParser(description='''
Produces the matrix branching program for conjunction on a given bitstring,
where 'bitstring' must contain only characters 0, 1, and ?.
Example: conjunction.py "01??1?11101?01"''')
    parser.add_argument('--cryfsm', metavar='PATH', action='store', type=str,
                        default='cryfsm',
                        help='set PATH as cryfsm executable')
    parser.add_argument('bitstring', action='store', type=str,
                        help='the bitstring to obfuscate')

    args = parser.parse_args()
    bitstring = args.bitstring
    r = re.search("[^01?]", bitstring)
    if r:
        print("error: invalid character '%c' in bitstring" % r.group())
        sys.exit(1)
    length = len(bitstring)
    fname = '%s.json' % re.escape(bitstring).replace('\\', '')
    lst = [args.cryfsm, 'conjunction.cry', '-g', '#', '-e', 'main \"%s\"' % bitstring,
           '-o', fname]
    try:
        run(lst)
    except OSError:
        print("error: failure when executing cryfsm")
        sys.exit(1)
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
