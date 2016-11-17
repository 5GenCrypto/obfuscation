#!/usr/bin/env python2

from __future__ import print_function
import os, re, subprocess, sys

def usage(ret):
    print("""
Usage: conjunction.py <bitstring>
    
Produces the matrix branching program for conjunction on a given bitstring,
where 'bitstring' must contain only characters 0, 1, and ?.
Example: conjunction.py 01??1?11101?01""")
    sys.exit(ret)

def run(lst):
    print('%s' % ' '.join(lst))
    with open(os.devnull, 'w') as fnull:
        return subprocess.call(lst, stdout=fnull, stderr=fnull)

def main(argv):
    if len(argv) != 2:
        usage(1)
    bitstring = argv[1]
    r = re.search("[^01?]", bitstring)
    if r:
        print("Error: invalid character '%c' in bitstring" % r.group())
        usage(1)
    fname = '%s.json' % re.escape(bitstring).replace('\\', '')
    lst = ['cryfsm', 'conjunction.cry', '-g', '#', '-e', 'main \"%s\"' % bitstring,
           '-o', fname]
    run(lst)

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
