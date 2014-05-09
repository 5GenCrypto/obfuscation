#!/usr/bin/env python2

from __future__ import division, print_function
import os, sys

def extract(line):
    return float(line.rstrip().rsplit(' ', 1)[1])

def main(argv):
    if len(argv) != 2:
        print('Usage: %s <path>' % argv[0])
        sys.exit(1)

    path = argv[1]

    obf_count = 0
    obf_num = 0

    mlm = False

    with open(path) as f:
        while True:
            try:
                line = f.next().strip()
            except StopIteration:
                break
            if line.startswith('Generating MLM'):
                mlm = True
            elif line.startswith('Generating p_i'):
                print('p_i/g_i: %f' % round(extract(line), 1))
            elif line.startswith('Generating CRT'):
                print('CRT: %f' % round(extract(line), 1))
            elif line.startswith('Generating z_i'):
                print('z_i: %f' % round(extract(line), 1))
            elif line.startswith('Generating pzt'):
                print('pzt: %f' % round(extract(line), 1))
            elif line.startswith('Took') and mlm:
                print('MLM Total: %f' % round(extract(line), 1))
                mlm = False
            elif line.startswith('Randomizing BPs'):
                line = f.next().strip()
                print('BP rand: %f' % round(extract(line), 1))
            elif line.startswith('Constructing bookend'):
                _ = f.next()
                a = f.next()
                _ = f.next()
                b = f.next()
                _ = f.next()
                total = extract(f.next())
                a = extract(a)
                b = extract(b)
                avg = (a + b) / 2
                print('Bookend enc (avg): %f' % round(avg, 1))
                print('Bookend enc (total): %f' % round(total, 1))
            elif line.startswith('Obfuscating layer'):
                _ = f.next()
                _ = f.next()
                a = f.next().strip()
                obf_count += extract(a)
                obf_num += 1
            elif line.startswith('Obfuscation took'):
                print('Layer enc (avg): %f' % round(obf_count / obf_num, 1))
                print('Layer enc (total): %f' % round(obf_count, 1))
                print('Obf Total: %f' % round(extract(line), 1))
            elif line.startswith('Max memory'):
                print('RAM (KB): %f' % round(extract(line) / 1024 / 1024, 1))
                break
                

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
