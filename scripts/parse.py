#!/usr/bin/env python2

import argparse, os, sys

def extract_time(line):
    return float(line.rsplit(' ', 1)[-1])
    
def parse_time(fname):
    with open(fname, 'r') as f:
        while True:
            try:
                line = f.next()
            except StopIteration:
                break
            if line.startswith('Generating MLM'):
                seed = f.next()
                primes = extract_time(f.next())
                x0 = extract_time(f.next())
                g_conversion = extract_time(f.next())
                crt = extract_time(f.next())
                zs = extract_time(f.next())
                pzt = extract_time(f.next())
                memoryusage = extract_time(f.next())
                saving = extract_time(f.next())
                total = extract_time(f.next())
            if line.startswith('Randomizing BPs'):
                total = extract_time(f.next())
                # Bookend vector stuff
                bookends = extract_time(f.next())
                first_encoding = extract_time(f.next())
                first_saving = extract_time(f.next())
                second_encoding = extract_time(f.next())
                second_saving = extract_time(f.next())
                total = extract_time(f.next())

def main(argv):
    parser = argparse.ArgumentParser(
        description='Parser for experiments.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--directory', action='store', default='runs',
                        help='directory containing run logs')

    args = parser.parse_args()
    for dir in os.listdir(args.directory):
        dirpath = os.path.join(args.directory, dir)
        for file in os.listdir(dirpath):
            path = os.path.join(dirpath, file)
            print(path)
            if 'time' in file:
                parse_time(path)


if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
