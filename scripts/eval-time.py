#!/usr/bin/env python2

import os, sys

import pylab
import numpy as np
import utils

def main(argv):
    utils.init()

    xs, idtimes = utils.dir_evaltime(os.path.join('results', 'secparam.id.circ'))
    _, andtimes = utils.dir_evaltime(os.path.join('results', 'secparam.and.circ'))
    _, xortimes = utils.dir_evaltime(os.path.join('results', 'secparam.xor.circ'))

    xs, idtimes, andtimes, xortimes = xs[1:], idtimes[1:], andtimes[1:], xortimes[1:]

    idtimes = [i * 1000 for i in idtimes]
    andtimes = [i * 1000 for i in andtimes]
    xortimes = [i * 1000 for i in xortimes]

    ind = np.arange(len(xs))
    width = 0.2

    pylab.figure(1)
    pylab.clf()
    pylab.axes([0.125,0.2,0.95-0.125,0.95-0.2])
    p1 = pylab.bar(ind + width, idtimes, width, color='black', log=True)
    p2 = pylab.bar(ind + 2 * width, andtimes, width, color='gray', log=True)
    p3 = pylab.bar(ind + 3 * width, xortimes, width, color='white', log=True)
    pylab.legend((p1[0], p2[0], p3[0]), ('ID', 'AND', 'XOR'), loc='upper left')
    pylab.xlabel(r'Security Parameter')
    pylab.ylabel(r'Time (ms)')
    pylab.xticks(np.arange(len(xs)) + 0.5, xs)
    pylab.ylim(100, 10 ** 4.2)

    pylab.savefig('/home/amaloz/umd-svn/iO-implementation/writeup/figs/single-gate-eval-time.eps')
    pylab.show()

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
