#!/usr/bin/env python2

from __future__ import division, print_function
import os, sys
import utils

def main(argv):
    if len(argv) != 2:
        print('Usage: %s <path>' % argv[0])
        sys.exit(1)

    path = argv[1]

    _, size = utils.obfsize(path)
    print('Size (MB): %f' % round(size / 2 ** 20, 5))
                

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
