#!/usr/bin/env python2

from __future__ import division

import os, sys

import pylab
import numpy as np
import utils

def extract(line):
    return float(line.rstrip().rsplit(' ', 1)[1])

def extract_obf_time(path):
    lst = []
    with open(path) as f:
        obf_count = 0
        obf_num = 0
        while True:
            try:
                line = f.next().strip()
            except StopIteration:
                break
            if line.startswith('Generating MLM'):
                mlm = True
            elif line.startswith('Generating p_i'):
                lst.append(round(extract(line), 1))
            elif line.startswith('Generating CRT'):
                lst.append(round(extract(line), 1))
            elif line.startswith('Generating z_i'):
                lst.append(round(extract(line), 1))
            elif line.startswith('Generating pzt'):
                lst.append(round(extract(line), 1))
            elif line.startswith('Took') and mlm:
                lst.append(round(extract(line), 1))
                mlm = False
            elif line.startswith('Randomizing BPs'):
                line = f.next().strip()
                lst.append(round(extract(line), 1))
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
                # lst.append(round(avg, 1))
                lst.append(round(total, 1))
            elif line.startswith('Obfuscating layer'):
                _ = f.next()
                _ = f.next()
                a = f.next().strip()
                obf_count += extract(a)
                obf_num += 1
            elif line.startswith('Obfuscation took'):
                # lst.append(round(obf_count / obf_num, 1))
                lst.append(round(obf_count, 1))
                lst.append(round(extract(line), 1))
            elif line.startswith('Max memory'):
                lst.append(round(extract(line) / 1024 / 1024, 1))
                break
    return lst

def extract_obf_size(path):
    _, size = utils.obfsize(path)
    return size

def extract_eval_time(path):
    _, time = utils.evaltime(path)
    return time

def compare_obf_time(point):
    a = extract_obf_time(os.path.join('results', 'point.52',
                                      'point-%d.circ-52-obf-time.log' % point))
    b = extract_obf_time(os.path.join('results', 'point.52.16cores',
                                      'point-%d.circ-52-obf-time.log' % point))
    p = [round((1 - i / j) * 100, 1) for i, j in zip(a, b)]
    print('Obf time: %d: %s' % (point, p))

def compare_obf_size(point):
    a = extract_obf_size(os.path.join('results', 'point.52',
                                      'point-%d.circ-52-obf-size.log' % point))
    b = extract_obf_size(os.path.join('results', 'point.52.16cores',
                                      'point-%d.circ-52-obf-size.log' % point))
    p = round((1 - a / b) * 100, 1)
    print('Obf size: %d: %s' % (point, p))

def compare_eval_time(point):
    a = extract_eval_time(os.path.join('results', 'point.52',
                                       'point-%d.circ-52-eval-time.log' % point))
    b = extract_eval_time(os.path.join('results', 'point.52.16cores',
                                       'point-%d.circ-52-eval-time.log' % point))
    p = round((1 - a / b) * 100, 1)
    print('Eval time: %d: %s' % (point, p))


def main(argv):
    compare_obf_time(8)
    compare_obf_size(8)
    compare_eval_time(8)
    compare_obf_time(12)
    compare_obf_size(12)
    compare_eval_time(12)
    compare_obf_time(16)
    compare_obf_size(16)
    compare_eval_time(16)

if __name__ == '__main__':
    try:
        main(sys.argv)
    except KeyboardInterrupt:
        pass
