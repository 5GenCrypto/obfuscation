#!/usr/bin/env python2

import os, sys

import matplotlib.pyplot as plt
import numpy as np
import utils

MB = 2 ** 20

def main(argv):
    utils.init_wide()

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

    _, idsizes = utils.dir_obfsize(idfiles)
    _, andsizes = utils.dir_obfsize(andfiles)
    _, xorsizes = utils.dir_obfsize(xorfiles)
    idsizes, andsizes, xorsizes = idsizes[1:], andsizes[1:], xorsizes[1:]
    idsizes = [i / MB for i in idsizes]
    andsizes = [i / MB for i in andsizes]
    xorsizes = [i / MB for i in xorsizes]

    _, ideval = utils.dir_evaltime(idfiles)
    _, andeval = utils.dir_evaltime(andfiles)
    _, xoreval = utils.dir_evaltime(xorfiles)
    ideval, andeval, xoreval = ideval[1:], andeval[1:], xoreval[1:]

    sizeexp = [(4 ** 5.39) * (x ** 2) for x in xs]
    evalexp = [(4 ** 5.84) * (x ** 2) / 1000 for x in xs]

    ind = np.arange(len(xs))
    width = 0.2

    fig, axes = plt.subplots(nrows=1, ncols=3)

    ax1 = axes.flat[0]

    p1 = ax1.bar(ind + width, idtimes, width, color='black')
    p2 = ax1.bar(ind + 2 * width, andtimes, width, color='gray')
    p3 = ax1.bar(ind + 3 * width, xortimes, width, color='white')
    ax1.legend((p1[0], p2[0], p3[0]), ('ID', 'AND', 'XOR'), loc='upper left')
    ax1.set_xlabel(r'Security Parameter')
    ax1.set_ylabel(r'Time (minutes)')
    ax1.set_xticks(ind + 0.5)
    ax1.set_xticklabels(xs)

    ax1 = axes.flat[1]

    p1 = ax1.bar(ind + width, idsizes, width, color='black')
    p2 = ax1.bar(ind + 2 * width, andsizes, width, color='gray')
    p3 = ax1.bar(ind + 3 * width, xorsizes, width, color='white')

    ax1.set_xlabel(r'Security Parameter')
    ax1.set_ylabel(r'Size (MB)')
    ax1.set_xticks(ind + 0.5)
    ax1.set_xticklabels(xs)

    ax2 = ax1.twinx()
    ax2.get_yaxis().set_visible(False)
    p4 = ax2.plot(ind + 2 * width + width / 2, sizeexp, 'k--')

    ax1.legend((p1[0], p2[0], p3[0]),
               ('ID', 'AND', 'XOR'),
               loc='upper left')

    ax1 = axes.flat[2]

    p1 = ax1.bar(ind + width, ideval, width, color='black')
    p2 = ax1.bar(ind + 2 * width, andeval, width, color='gray')
    p3 = ax1.bar(ind + 3 * width, xoreval, width, color='white')


    ax1.set_ylabel(r'Time (s)')
    ax1.set_xlabel(r'Security Parameter')
    ax1.set_xticks(ind + 0.5)
    ax1.set_xticklabels(xs)

    ax2 = ax1.twinx()
    ax2.get_yaxis().set_visible(False)
    p4 = ax2.plot(ind + 2 * width + width / 2, evalexp, 'k--')

    ax1.legend((p1[0], p2[0], p3[0]),
               ('ID', 'AND', 'XOR'),
               loc='upper left')

    plt.subplots_adjust(wspace=0.4, bottom=0.2)
    plt.savefig('/home/amaloz/umd-svn/iO-implementation/writeup/figs/secparam.eps')
    plt.show()

    # pylab.figure(1)
    # pylab.clf()
    # pylab.axes([0.125,0.2,0.95-0.125,0.95-0.2])
    # p1 = pylab.bar(ind + width, idtimes, width, color='black')
    # p2 = pylab.bar(ind + 2 * width, andtimes, width, color='gray')
    # p3 = pylab.bar(ind + 3 * width, xortimes, width, color='white')
    # pylab.legend((p1[0], p2[0], p3[0]), ('ID', 'AND', 'XOR'), loc='upper left')
    # pylab.xlabel(r'Security Parameter')
    # pylab.ylabel(r'Time (minutes)')
    # pylab.xticks(np.arange(len(xs)) + 0.5, xs)

    # pylab.savefig('/home/amaloz/umd-svn/iO-implementation/writeup/figs/single-gate-obf-time.eps')
    # pylab.show()

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
