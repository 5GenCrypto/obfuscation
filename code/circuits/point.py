#!/usr/bin/env python2

import random, sys

def error():
    print('Usage: point.py <bitlength>')
    sys.exit(1)

def random_bitstring(bitlength):
    r = random.randint(0, 2 ** bitlength - 1)
    bits = bin(r)[2:]
    return '0' * (bitlength - len(bits)) + bits

def main(argv):
    if len(argv) != 2:
        error()
    try:
        bitlength = int(argv[1])
    except ValueError:
        error()

    with open('point-%d.circ' % bitlength, 'w') as f:
        secret = random_bitstring(bitlength)
        for _ in xrange(5):
            test = random_bitstring(bitlength)
            if test != secret:
                f.write('# TEST %s 0\n' % test)
        f.write('# TEST %s 1\n' % secret)
        for i in xrange(bitlength):
            f.write('%d input\n' % i)
        for i in xrange(bitlength):
            gate = 'ID' if secret[i] == '0' else 'NOT'
            f.write('%d gate %s %d\n' % (bitlength + i, gate, i))

        length = bitlength
        oldstart = bitlength
        start = oldstart + length
        while length > 1:
            for i, j in zip(xrange(start, start + length / 2),
                            xrange(oldstart, start, 2)):
                f.write('%d gate OR %d %d\n' % (i, j, j + 1))
            length = length / 2
            oldstart = start
            start = start + length
        f.write('%d output NOT %d\n' % (start, start - 1))

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
