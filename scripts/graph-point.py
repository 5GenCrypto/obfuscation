#!/usr/bin/env python2

from __future__ import division
import os, sys

from operator import add, sub
import matplotlib.pyplot as plt
import numpy as np
import utils

def obftime(ax1):
    xs = (8, 12, 16)

    path8 = os.path.join('results', 'point.52', 'point-8.circ-52-obf-time.log')
    path12 = os.path.join('results', 'point.52', 'point-12.circ-52-obf-time.log')
    path16 = os.path.join('results', 'point.52', 'point-16.circ-52-obf-time.log')

    bar8 = utils.extract_obf_time(path8)
    bar12 = utils.extract_obf_time(path12)
    bar16 = utils.extract_obf_time(path16)

    all = zip(bar8, bar12, bar16)
    all = [(a / 60 / 60, b / 60 / 60, c / 60 / 60) for a, b, c in all]

    ind = np.arange(len(all[0]))
    width = 0.4

    expected = [(l ** 2) * (s ** 5) * (52 ** 2)
                for l, s in zip(xs, map(sub, xs, [1, 1, 1]))]

    colors = [(0, i / len(all), 0, 1) for i in xrange(len(all))]

    encodingtime = map(add, all[6], all[7])
    # other = map(sub, all[8], encodingtime)
    # other = map(sub, other, all[4])

    # ax1 = plt.axes([0.125,0.2,0.95-0.125,0.95-0.2])
    total = ax1.bar(ind + width, all[8], width, color='white')
    mlm = ax1.bar(ind + width, all[4], width, color='black')
    enc = ax1.bar(ind + width, encodingtime, width, color='gray', bottom=all[4])
    # pigi = ax1.bar(ind + width, all[0], width, color=colors[0])
    # psum = all[0]
    # crt = ax1.bar(ind + width, all[1], width, color=colors[1], bottom=psum)
    # psum = map(add, all[1], psum)
    # zi = ax1.bar(ind + width, all[2], width, color=colors[2], bottom=psum)
    # psum = map(add, all[2], psum)
    # pzt = ax1.bar(ind + width, all[3], width, color=colors[3], bottom=psum)
    # psum = all[4]
    # bprand = ax1.bar(ind + width, all[5], width, color=colors[5], bottom=psum)
    # psum = map(add, all[5], psum)
    # bookends = ax1.bar(ind + width, all[6], width, color=colors[6], bottom=psum)
    # psum = map(add, all[6], psum)
    # layers = ax1.bar(ind + width, all[7], width, color=colors[7], bottom=psum)

    ax1.set_ylabel(r'Obfuscation time (hr)')
    ax1.set_xlabel(r'Input size of point function')
    ax1.set_xticks(ind + 0.6)
    ax1.set_ylim(0, 10)
    ax1.set_xticklabels(xs)

    # ax2 = ax1.twinx()
    # ax2.get_yaxis().set_visible(False)
    # expline = ax2.plot(ind + 0.6, expected, 'k--')

    ax1.legend((mlm[0], enc[0]),
               ('Param gen', 'Encoding'),
               loc='upper left')

def obfsize(ax1):
    xs = (8, 12, 16)

    path8 = os.path.join('results', 'point.52', 'point-8.circ-52-obf-size.log')
    path12 = os.path.join('results', 'point.52', 'point-12.circ-52-obf-size.log')
    path16 = os.path.join('results', 'point.52', 'point-16.circ-52-obf-size.log')

    _, bar8 = utils.obfsize(path8)
    _, bar12 = utils.obfsize(path12)
    _, bar16 = utils.obfsize(path16)

    all = [bar8, bar12, bar16]
    all = [x / 2 ** 30 for x in all]

    ind = np.arange(len(all))
    width = 0.4

    expected = [(l ** 2) * (s ** 5) * (52**2)
                for l, s in zip(xs, map(sub, xs, [1, 1, 1]))]

    total = ax1.bar(ind + width, all, width, color='gray')
    ax1.set_ylabel(r'Obfuscation size (GB)')
    ax1.set_xlabel(r'Input size of point function')
    ax1.set_xticks(ind + 0.6)
    ax1.set_xticklabels(xs)

    ax2 = ax1.twinx()
    ax2.get_yaxis().set_visible(False)
    expline = ax2.plot(ind + 0.6, expected, 'k--')

    # ax1.legend((expline[0],),
    #            ('$O(n^2 s^5 \lambda^2)$',),
    #            loc='upper left')

def evaltime(ax1):
    xs = (8, 12, 16)

    path8 = os.path.join('results', 'point.52', 'point-8.circ-52-eval-time.log')
    path12 = os.path.join('results', 'point.52', 'point-12.circ-52-eval-time.log')
    path16 = os.path.join('results', 'point.52', 'point-16.circ-52-eval-time.log')

    _, bar8 = utils.evaltime(path8)
    _, bar12 = utils.evaltime(path12)
    _, bar16 = utils.evaltime(path16)

    all = [bar8, bar12, bar16]
    all = [x / 60 / 60 for x in all]

    ind = np.arange(len(all))
    width = 0.4

    expected = [(l ** 2) * (s ** 5.37) * (52**2)
                for l, s in zip(xs, map(sub, xs, [1, 1, 1]))]

    total = ax1.bar(ind + width, all, width, color='gray')
    ax1.set_ylabel(r'Evaluation time (hr)')
    ax1.set_xlabel(r'Input size of point function')
    ax1.set_xticks(ind + 0.6)
    ax1.set_xticklabels(xs)

    ax2 = ax1.twinx()
    ax2.get_yaxis().set_visible(False)
    expline = ax2.plot(ind + 0.6, expected, 'k--')

    # ax1.legend((expline[0],),
    #            ('$O(n^2 s^{5.37} \lambda^2)$',),
    #            loc='upper left')


    # plt.savefig('/home/amaloz/umd-svn/iO-implementation/writeup/figs/point-eval-time.eps')
    # plt.show()


def main(argv):
    utils.init_wide()

    fig, axes = plt.subplots(nrows=1, ncols=3)
    obftime(axes.flat[0])
    obfsize(axes.flat[1])
    evaltime(axes.flat[2])

    plt.subplots_adjust(wspace=0.4, bottom=0.2)
    plt.savefig('/home/amaloz/umd-svn/iO-implementation/writeup/figs/point-obf.eps')
    plt.show()

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
