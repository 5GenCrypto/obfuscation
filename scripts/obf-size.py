#!/usr/bin/env python2

from __future__ import division
import os, sys

import matplotlib.pyplot as plt
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

    expected = [(4 ** 5.39) * (x ** 2) for x in xs]
    print(expected)

    ax1 = plt.axes([0.125,0.2,0.95-0.125,0.95-0.2])
    p1 = ax1.bar(ind + width, idsizes, width, color='black')
    p2 = ax1.bar(ind + 2 * width, andsizes, width, color='gray')
    p3 = ax1.bar(ind + 3 * width, xorsizes, width, color='white')

    ax1.set_xlabel(r'Security Parameter')
    ax1.set_ylabel(r'Size (MB)')
    ax1.set_xticks(ind + 0.5)
    ax1.set_xticklabels(xs)

    ax2 = ax1.twinx()
    ax2.get_yaxis().set_visible(False)
    p4 = ax2.plot(ind + 2 * width + width / 2, expected, 'k--')

    ax1.legend((p1[0], p2[0], p3[0], p4[0]),
               ('ID', 'AND', 'XOR', '$\mathcal{O}(4^{5.39}\lambda^2)$'),
               loc='upper left')
    plt.savefig('/home/amaloz/umd-svn/iO-implementation/writeup/figs/single-gate-obf-size.eps')
    plt.show()

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
