#!/usr/bin/env python2

from __future__ import division

import os, sys

import pylab
import numpy as np
import utils

MB = 2 ** 20

def main(argv):
    utils.init()

    xs, idsizes = utils.dir_obfsize(os.path.join('results', 'secparam.id.circ'))
    _, andsizes = utils.dir_obfsize(os.path.join('results', 'secparam.and.circ'))
    _, xorsizes = utils.dir_obfsize(os.path.join('results', 'secparam.xor.circ'))

    xs, idsizes, andsizes, xorsizes = xs[1:], idsizes[1:], andsizes[1:], xorsizes[1:]

    idsizes = [i / MB for i in idsizes]
    andsizes = [i / MB for i in andsizes]
    xorsizes = [i / MB for i in xorsizes]

    print(idsizes)
    print(andsizes)
    print(xorsizes)

    ind = np.arange(len(xs))
    width = 0.2

    pylab.figure(1)
    pylab.clf()
    pylab.axes([0.125,0.2,0.95-0.125,0.95-0.2])
    p1 = pylab.bar(ind + width, idsizes, width, color='black', log=True)
    p2 = pylab.bar(ind + 2 * width, andsizes, width, color='gray', log=True)
    p3 = pylab.bar(ind + 3 * width, xorsizes, width, color='white', log=True)
    pylab.legend((p1[0], p2[0], p3[0]), ('ID', 'AND', 'XOR'), loc='upper left')
    pylab.xlabel(r'Security Parameter')
    pylab.ylabel(r'Size (MB)')
    pylab.xticks(np.arange(len(xs)) + 0.5, xs)
    pylab.ylim(1, 10 ** 3)

    pylab.savefig('/home/amaloz/umd-svn/iO-implementation/writeup/figs/single-gate-obf-size.eps')
    pylab.show()

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
