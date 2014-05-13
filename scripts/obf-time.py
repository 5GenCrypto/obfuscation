#!/usr/bin/env python2

import os, sys

import pylab
import numpy as np
import utils

def main(argv):
    utils.init()

    idfiles = os.path.join('results', 'secparam.id.circ')
    andfiles = os.path.join('results', 'secparam.and.circ')
    xorfiles = os.path.join('results', 'secparam.xor.circ')
    
    xs, idtimes = utils.dir_obftime(idfiles)
    _, andtimes = utils.dir_obftime(andfiles)
    _, xortimes = utils.dir_obftime(xorfiles)

    xs, idtimes, andtimes, xortimes = xs[1:], idtimes[1:], andtimes[1:], xortimes[1:]

    idtimes = [x / 60 for x in idtimes]
    andtimes = [x / 60 for x in andtimes]
    xortimes = [x / 60 for x in xortimes]

    ind = np.arange(len(xs))
    width = 0.2

    pylab.figure(1)
    pylab.clf()
    pylab.axes([0.125,0.2,0.95-0.125,0.95-0.2])
    p1 = pylab.bar(ind + width, idtimes, width, color='black')
    p2 = pylab.bar(ind + 2 * width, andtimes, width, color='gray')
    p3 = pylab.bar(ind + 3 * width, xortimes, width, color='white')
    pylab.legend((p1[0], p2[0], p3[0]), ('ID', 'AND', 'XOR'), loc='upper left')
    pylab.xlabel(r'Security Parameter')
    pylab.ylabel(r'Time (minutes)')
    pylab.xticks(np.arange(len(xs)) + 0.5, xs)
    # pylab.ylim(1, 10 ** 3.2)

    pylab.savefig('/home/amaloz/umd-svn/iO-implementation/writeup/figs/single-gate-obf-time.eps')
    pylab.show()

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
