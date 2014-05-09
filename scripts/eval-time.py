#!/usr/bin/env python2

import os, sys

import pylab
import numpy as np
import utils

def main(argv):
    utils.init()

    xs, id_y = utils.dir_evaltime(os.path.join('results', 'secparam.id.circ'))
    _, and_y = utils.dir_evaltime(os.path.join('results', 'secparam.and.circ'))
    _, xor_y = utils.dir_evaltime(os.path.join('results', 'secparam.xor.circ'))

    xs = xs[::2]
    id_y = id_y[::2]
    and_y = and_y[::2]
    xor_y = xor_y[::2]

    id_y = [i * 1000.0 for i in id_y]
    and_y = [i * 1000.0 for i in and_y]
    xor_y = [i * 1000.0 for i in xor_y]

    ind = np.arange(len(xs))
    width = 0.2

    pylab.figure(1)
    pylab.clf()
    pylab.axes([0.125,0.2,0.95-0.125,0.95-0.2])
    p1 = pylab.bar(ind + width, id_y, width, color='black', log=True)
    p2 = pylab.bar(ind + 2 * width, and_y, width, color='gray', log=True)
    p3 = pylab.bar(ind + 3 * width, xor_y, width, color='white', log=True)
    pylab.legend((p1[0], p2[0], p3[0]), ('ID', 'AND', 'XOR'), loc='upper left')
    pylab.xlabel(r'Security Parameter')
    pylab.ylabel(r'Time (ms)')
    pylab.xticks(np.arange(len(xs)) + 0.5, xs)
    # pylab.ylim(1, 10 ** 3.2)

    pylab.savefig('/home/amaloz/umd-svn/iO-implementation/writeup/figs/single-gate-eval-time.eps')
    pylab.show()

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
